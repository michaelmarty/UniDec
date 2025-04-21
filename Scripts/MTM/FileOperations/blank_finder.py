import os
import shutil
def find_blanks(topdir):
    blanks = []
    for root, dirs, files in os.walk(topdir):
        for d in dirs:
            if "blank" in d.lower() or "wash" in d.lower():
                blanks.append(os.path.join(root, d))
    return blanks

def delete_blanks(blanks):
    for b in blanks:
        print("Deleting:", b)
        # rmtree
        shutil.rmtree(b)


topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Melanie Odenkirk"
print("Searching for blanks in:", topdir)
blanks = find_blanks(topdir)
print("Blanks Found:")
for b in blanks:
    print(b)

# delete_blanks(blanks)