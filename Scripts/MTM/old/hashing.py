import sys
import hashlib


def hashfile(path):
    # BUF_SIZE is totally arbitrary, change for your app!
    BUF_SIZE = 65536  # lets read stuff in 64kb chunks!

    md5 = hashlib.md5()
    sha256 = hashlib.sha256()

    with open(path, 'rb') as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            md5.update(data)
            sha256.update(data)

    md5hash = md5.hexdigest()
    sha256hash = sha256.hexdigest()

    print("MD5: {0}".format(md5hash))
    print("SHA256: {0}".format(sha256hash))

    return md5hash, sha256hash

if __name__ == "__main__":
    hashfile("dist\\UniDec_Windows_230920.zip")