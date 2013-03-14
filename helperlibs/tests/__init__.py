from os import path
_TESTFILE_DIR='files'
def get_file_path(filename):
    "Get the path of the test file"
    return path.join(path.abspath(path.dirname(__file__)),_TESTFILE_DIR, filename)
