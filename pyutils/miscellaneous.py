from pyutils import DATA_DIR, logging
from subprocess import Popen, CalledProcessError, PIPE, STDOUT


def run_cmd(*args):
    executable = args[0]
    with Popen(
        args=args,
        stdout=PIPE,
        stderr=STDOUT,
        cwd=DATA_DIR,
        bufsize=1,
        universal_newlines=True
    ) as p:
        for line in p.stdout:
            logging.info(f"[{executable}]: {line.strip()}")
        if p.returncode:
            raise CalledProcessError(p.returncode, p.args)
