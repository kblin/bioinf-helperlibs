import tempfile
import shutil
import os
from os import path
import sys

class TemporaryFile(object):
    def __init__(self, suffix='', prefix='tmp', dir=None, text=False):
        self.handle, self.name = tempfile.mkstemp(suffix, prefix, dir, text)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if hasattr(self.handle, "close"):
            self.handle.close()
        os.unlink(self.name)


class TemporaryDirectory(object):
    def __init__(self, suffix='', prefix='tmp', dir=None, change=False):
        self.change = change
        self.tempdir = tempfile.mkdtemp(suffix, prefix, dir)

    def __enter__(self):
        if self.change:
            self.old_wd = os.getcwd()
            os.chdir(self.tempdir)
        return self.tempdir

    def __exit__(self, type, value, traceback):
        if self.change:
            os.chdir(self.old_wd)
        shutil.rmtree(self.tempdir)


if sys.platform in ('linux2', 'darwin'):
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
else:
    class TemporaryPipe(object):
        def __init__(self, pipename=''):
            self.handle, self.name = tempfile.mkstemp()

        def __enter__(self):
            return self.name

        def __exit__(self, type, value, traceback):
            if hasattr(self.handle, "close"):
                self.handle.close()
            os.unlink(self.name)
