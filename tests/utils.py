import os

def exists(p, msg):
    assert os.path.exists(p), msg

def list_files(in_path, ext = None):
    files = []
    for (dirpath, dirnames, filenames) in os.walk(in_path):
        if ext is not None:
            files.extend([f for f in filenames if f.rfind(ext) > 0])
        else:
            files.extend(filenames)
        break

    return files

def get_files_from_dir(path, ext = None):
    files = list_files(path, ext = ext)
    return [os.path.join(path, x) for x in files]