"""Regression tests for the release workflow; no Magma or network access needed."""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
TIME_LINE = re.compile(r'^\s*[0-9]+\.[0-9]+\s*$')


def closing_line(text):
    """The closing line of a log, with the rules of verify_all.sh."""
    lines = [line.rstrip('\r') for line in text.splitlines()]
    lines = [line for line in lines if line.strip()]
    while lines and (lines[-1].startswith('Script done on ') or TIME_LINE.match(lines[-1])):
        lines.pop()
    return lines[-1] if lines else ''


def manifest():
    rows = []
    for line in (ROOT / 'jobs.txt').read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        rows.append(line.split('\t'))
    return rows


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        shutil.copy2(ROOT / 'verify_all.sh', self.root)
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        self.env = dict(os.environ, PATH=str(self.bin) + os.pathsep + os.environ['PATH'])

    def jobs(self, *rows):
        (self.root / 'jobs.txt').write_text('\n'.join('\t'.join(row) for row in rows) + '\n')

    def run_jobs(self, *args):
        return subprocess.run(['bash', './verify_all.sh', *args], cwd=self.root, env=self.env,
                              capture_output=True, text=True, timeout=30)

    def fake_magma(self):
        for name in ('mdmagma/v2/mdmagma.spec', 'mdmagma/Magma/magma.spec',
                     'gl2/magma.spec', 'fast_hecke.m', 'mdmagma-magma2.30.patch'):
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        fake = self.bin / 'magma'
        fake.write_text('''#!/usr/bin/env python3
import json, os, pathlib, sys
with open('magma_calls.jsonl', 'a') as f:
    f.write(json.dumps(sys.argv[1:]) + '\\n')
if pathlib.Path(sys.argv[-1]).name == 'init.m':
    print(os.environ.get('INIT_OUTPUT', 'PACKAGE_READY'))
    sys.exit(int(os.environ.get('INIT_EXIT', '0')))
print(os.environ.get('JOB_OUTPUT', 'DONE'))
''')
        fake.chmod(0o755)
        # git is only used to check and apply the mdmagma patch
        git = self.bin / 'git'
        git.write_text('''#!/usr/bin/env python3
import json, os, sys
with open('git_calls.jsonl', 'a') as f:
    f.write(json.dumps(sys.argv[1:]) + '\\n')
if '--reverse' in sys.argv:
    sys.exit(int(os.environ.get('GIT_REVERSE_CHECK', '0')))
sys.exit(int(os.environ.get('GIT_APPLY', '0')))
''')
        git.chmod(0o755)

    def calls(self, name):
        path = self.root / name
        if not path.exists():
            return []
        return [json.loads(line) for line in path.read_text().splitlines()]

    def test_serial_and_parallel_success(self):
        self.jobs(('logs/a.log', "printf 'RESULT survivors=0\\nDONE\\n\\n'", 'DONE',
                   '^RESULT survivors=0$'), ('logs/b.log', "printf 'DONE details\\n'", 'DONE'))
        for args in ((), ('-j', '2')):
            with self.subTest(args=args):
                result = self.run_jobs(*args)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn('all 2 jobs completed', result.stdout)

    def test_closing_line_ignores_cpu_times_carriage_returns_and_script_trailer(self):
        self.jobs(('logs/a.log', "printf 'SIEVE_DONE 48 0 []\\n13.890\\n'", 'SIEVE_DONE 48 0 []'),
                  ('logs/b.log', "printf 'DONE 1.5\\r\\n\\r\\nScript done on today\\n'", 'DONE'),
                  ('logs/c.log', "printf 'M3_CHECK 3\\n'", 'M3_CHECK 3'))
        result = self.run_jobs()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.jobs(('logs/a.log', "printf 'DONE\\nSOMETHING ELSE\\n13.890\\n'", 'DONE'))
        self.assertEqual(self.run_jobs().returncode, 1)

    def test_failed_commands_incomplete_logs_and_wrong_results_fail(self):
        failures = [
            ("bash -c 'echo DONE; exit 7'", 'DONE'),
            ("false | printf 'DONE\\n'", 'DONE'),
            ("exit 0", 'DONE'),
            ("printf 'Magma: Internal error\\n'", 'DONE'),
            ("printf 'DONE_WRONG\\n'", 'DONE'),
            ("printf 'RESULT survivors=3\\nDONE\\n'", 'DONE', '^RESULT survivors=0$'),
        ]
        for failure in failures:
            for args in ((), ('-j', '2')):
                with self.subTest(failure=failure, args=args):
                    self.jobs(('logs/failed.log', *failure),
                              ('logs/good.log', "printf 'DONE\\n'", 'DONE'))
                    result = self.run_jobs(*args)
                    self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                    self.assertIn('1 of 2 jobs failed', result.stdout)

    def test_job_failure_count(self):
        self.jobs(('logs/a.log', 'false', 'DONE'), ('logs/b.log', 'false', 'DONE'))
        self.assertEqual(self.run_jobs('-j', '2').returncode, 2)

    def test_invalid_invocations_do_not_overwrite_logs(self):
        self.jobs(('logs/a.log', "printf 'DONE\\n'", 'DONE'))
        (self.root / 'logs').mkdir()
        (self.root / 'logs/a.log').write_text('reference\n')
        for args in (('-j', '0'), ('-j', '-2'), ('-j', 'bad'), ('-j',),
                     ('missing',), ('-j', '2', 'missing'), ('a', 'b')):
            with self.subTest(args=args):
                self.assertNotEqual(self.run_jobs(*args).returncode, 0)
                self.assertEqual((self.root / 'logs/a.log').read_text(), 'reference\n')

    def test_missing_empty_and_malformed_manifests_fail(self):
        self.assertNotEqual(self.run_jobs().returncode, 0)
        for content in ('# no jobs\n', 'logs/a.log\ttrue\n',
                        'logs/a.log\ttrue\tDONE\t[\n',
                        'logs/a.log\ttrue\tDONE\nlogs/a.log\ttrue\tDONE\n'):
            with self.subTest(content=content):
                (self.root / 'jobs.txt').write_text(content)
                self.assertNotEqual(self.run_jobs().returncode, 0)

    def test_archive_only_jobs_are_not_run(self):
        (self.root / 'logs').mkdir()
        (self.root / 'logs/old.log').write_text('archived output\nPROGRESS 1/2\n')
        self.jobs(('logs/old.log', 'archived: interrupted run of: printf DONE', 'PROGRESS 1/2'),
                  ('logs/new.log', "printf 'DONE\\n'", 'DONE'))
        result = self.run_jobs()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('skipping archive-only logs/old.log', result.stdout)
        self.assertEqual((self.root / 'logs/old.log').read_text(), 'archived output\nPROGRESS 1/2\n')
        result = self.run_jobs('old')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((self.root / 'logs/old.log').read_text(), 'archived output\nPROGRESS 1/2\n')

    def test_python_jobs_need_no_magma_and_run_after_magma_jobs(self):
        # Without magma on PATH, a selection of non-Magma jobs still runs.
        self.jobs(('logs/py.log', "printf 'DONE\\n'", 'DONE'))
        env = dict(self.env, PATH=str(self.bin) + os.pathsep + '/usr/bin' + os.pathsep + '/bin')
        result = subprocess.run(['bash', './verify_all.sh'], cwd=self.root, env=env,
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        # A non-Magma job listed first still sees the output of the Magma job.
        self.fake_magma()
        self.env['JOB_OUTPUT'] = 'MAGMA_DONE'
        self.jobs(('logs/reader.log', 'cat logs/magma.log', 'MAGMA_DONE'),
                  ('logs/magma.log', 'magma -n -b proof.m', 'MAGMA_DONE'))
        for args in ((), ('-j', '2')):
            with self.subTest(args=args):
                result = self.run_jobs(*args)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_scheduler_cannot_silently_skip_jobs(self):
        self.fake_magma()
        self.jobs(('logs/a.log', 'magma -n -b proof.m', 'DONE'))
        fake = self.bin / 'xargs'
        for status in (0, 7):
            with self.subTest(status=status):
                fake.write_text('#!/bin/sh\nexit %d\n' % status)
                fake.chmod(0o755)
                self.assertNotEqual(self.run_jobs('-j', '2').returncode, 0)

    def test_package_initialisation_checks_status_and_exact_marker(self):
        self.fake_magma()
        self.jobs(('logs/a.log', 'magma -n -b proof.m', 'DONE'))
        for output, status in [('PACKAGE_READY', '7'),
                               ('> print "PACKAGE_READY";\nMagma: Internal error', '0'),
                               ('PACKAGE_READY_BROKEN', '0')]:
            with self.subTest(output=output, status=status):
                self.env.update(INIT_OUTPUT=output, INIT_EXIT=status)
                self.assertNotEqual(self.run_jobs('-j', '2').returncode, 0)
                self.assertFalse((self.root / 'logs/a.log').exists())

    def test_magma_error_in_log_fails_the_job(self):
        self.fake_magma()
        self.jobs(('logs/a.log', "printf 'Runtime error: x\\nDONE\\n'", 'DONE'))
        result = self.run_jobs()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('MAGMA ERROR', result.stdout)

    def test_packages_initialised_sequentially_before_jobs(self):
        self.fake_magma()
        self.jobs(('logs/a.log', 'magma -n -b proof.m', 'DONE'))
        result = self.run_jobs('-j', '2')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        calls = self.calls('magma_calls.jsonl')
        self.assertEqual(len(calls), 4)
        self.assertTrue(all(call[:2] == ['-n', '-b'] for call in calls))
        self.assertTrue(all(Path(call[-1]).name == 'init.m' for call in calls[:3]))
        self.assertEqual(calls[-1][-1], 'proof.m')

    def test_mdmagma_patch_is_applied_once(self):
        self.fake_magma()
        self.jobs(('logs/a.log', 'magma -n -b proof.m', 'DONE'))
        # already applied: only the reverse check runs
        result = self.run_jobs()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual([c for c in self.calls('git_calls.jsonl') if '--reverse' not in c], [])
        # not yet applied: the patch is applied
        (self.root / 'git_calls.jsonl').unlink()
        self.env['GIT_REVERSE_CHECK'] = '1'
        result = self.run_jobs()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(['-C', 'mdmagma', 'apply', '../mdmagma-magma2.30.patch'],
                      self.calls('git_calls.jsonl'))
        # cannot be applied: nothing runs
        (self.root / 'logs/a.log').unlink()
        self.env['GIT_APPLY'] = '1'
        self.assertNotEqual(self.run_jobs().returncode, 0)
        self.assertFalse((self.root / 'logs/a.log').exists())


class ArchiveTests(unittest.TestCase):
    def test_manifest_covers_logs_and_matches_reference_results(self):
        paths = []
        for fields in manifest():
            self.assertIn(len(fields), (3, 4))
            self.assertTrue(all(fields), fields)
            path, command, marker = fields[:3]
            paths.append(path)
            log = (ROOT / path).read_text(errors='replace')
            last = closing_line(log)
            self.assertTrue(last == marker or last.startswith(marker + ' '), path)
            if len(fields) == 4:
                self.assertEqual(subprocess.run(['grep', '-aEq', '--', fields[3]], input=log,
                                               text=True).returncode, 0, path)
        self.assertEqual(len(paths), len(set(paths)))
        self.assertEqual(set(paths), {str(path.relative_to(ROOT))
                                     for folder in ('logs', 'independent/logs')
                                     for path in (ROOT / folder).rglob('*.log')})

    def test_commands_refer_to_files_in_the_repository(self):
        for path, command, *_ in manifest():
            if command.startswith('archived:'):
                self.assertGreater(len(command), len('archived: '), path)
                continue
            words = command.split()
            self.assertIn(words[0], ('magma', 'python3'), path)
            for word in words:
                value = word.split(':=')[-1]
                if value.endswith(('.m', '.py', '.txt', '.json')):
                    self.assertTrue((ROOT / value).is_file(), (path, value))

    def test_magma_dependencies_are_present(self):
        pattern = re.compile(r'(?:AttachSpec|Attach|Read)\(\s*"([^"]+)"|load\s+"([^"]+)"')
        for script in ROOT.glob('*.m'):
            for match in pattern.finditer(script.read_text()):
                target = match.group(1) or match.group(2)
                if target.startswith(('mdmagma/', 'gl2/')):
                    continue  # submodules; checked by verify_all.sh before running Magma
                self.assertTrue((ROOT / target).is_file(), (script.name, target))


if __name__ == '__main__':
    unittest.main()
