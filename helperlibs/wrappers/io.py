import tempfile
import shutil
import os
from os import path

class TemporaryPipe(object):
    def __init__(self, pipename="pipe"):
        self.pipename = pipename
        self.tempdir = None

    def __enter__(self):
        self.tempdir = tempfile.mkdtemp()
        pipe_path = path.join(self.tempdir, self.pipename)
        os.mkfifo(pipe_path)
        return pipe_path

    def __exit__(self, type, value, traceback):
        if self.tempdir is not None:
            shutil.rmtree(self.tempdir)

class TemporaryFile(object):
    def __init__(self, suffix='', prefix='tmp', dir=None, text=False):
        self.handle, self.name = tempfile.mkstemp(suffix, prefix, dir, text)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if hasattr(self.handle, "close"):
            self.handle.close()
        os.unlink(self.name)
