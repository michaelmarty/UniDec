"""
Rename files in a folder replacing 'AmgenhumanmAb' with 'AmgenrabbitmAb'.

Usage (PowerShell):
	python AmgenhumanmAb_AmgenrabbitmAb_rename.py "Z:\\Group Share\\BHT\\Temp\\Paper\\DMT\\Amgen Rabbit IgG 21 Injections Set 1"

Options:
	--dry-run    : Show what would be renamed without performing changes
	--pattern    : Pattern to replace (default: AmgenhumanmAb)
	--replace    : Replacement string (default: AmgenrabbitmAb)

The script performs a case-insensitive search and replaces any occurrence of the pattern
in filenames with the replacement string (the replacement uses the exact case you pass).
"""

import os
import sys
import argparse
import re


def rename_in_folder(folder, pattern='AmgenhumanmAb', replacement='AmgenrabbitmAb', dry_run=True):
	folder = os.path.abspath(folder)
	if not os.path.isdir(folder):
		raise NotADirectoryError(f"Folder does not exist: {folder}")

	files = os.listdir(folder)
	changed = []
	for fname in files:
		# perform case-insensitive search
		if re.search(re.escape(pattern), fname, flags=re.IGNORECASE):
			newfname = re.sub(re.escape(pattern), replacement, fname, flags=re.IGNORECASE)
			if newfname == fname:
				continue
			src = os.path.join(folder, fname)
			dst = os.path.join(folder, newfname)
			if os.path.exists(dst):
				print(f"Skipping rename because destination exists: {dst}")
				continue
			changed.append((src, dst))

	if not changed:
		print("No files to rename.")
		return changed

	# Print planned changes
	print("Planned renames:")
	for src, dst in changed:
		print(f"  {os.path.basename(src)}  ->  {os.path.basename(dst)}")

	if dry_run:
		print("Dry run: no files were changed. Rerun without --dry-run to perform renames.")
		return changed

	# Perform renames
	for src, dst in changed:
		try:
			os.rename(src, dst)
			print(f"Renamed: {os.path.basename(src)} -> {os.path.basename(dst)}")
		except Exception as e:
			print(f"Failed to rename {src} -> {dst}: {e}")

	return changed


if __name__ == '__main__':
	p = argparse.ArgumentParser(description='Rename AmgenhumanmAb -> AmgenrabbitmAb in filenames')
	p.add_argument('folder', nargs='?', default=None,
				   help='Folder to scan (you can paste a Windows path)')
	# If these are not provided on the command line the script will prompt for them
	p.add_argument('--pattern', default=None, help="Pattern to look for (case-insensitive). If omitted you'll be prompted")
	p.add_argument('--replace', dest='replace', default=None, help="Replacement string. If omitted you'll be prompted")
	# --apply when present means perform renames; otherwise the script will prompt
	p.add_argument('--apply', dest='apply', action='store_true', default=False, help='If set, perform the renames; otherwise you will be prompted')

	args = p.parse_args()

	# If no folder was provided on the command line, prompt the user to paste one.
	folder = args.folder
	if not folder:
		try:
			folder = input('Enter folder path (paste Windows path and press Enter): ').strip()
		except EOFError:
			print('No folder provided; exiting.')
			sys.exit(1)

	# Strip surrounding quotes and expand environment variables/user (~)
	folder = folder.strip('"').strip("'")
	folder = os.path.expanduser(os.path.expandvars(folder))

	# Prompt for options if they were not provided on the command line
	if args.pattern is None:
		pat = input("Pattern to replace [default: AmgenhumanmAb]: ").strip()
		pattern = pat or 'AmgenhumanmAb'
	else:
		pattern = args.pattern

	if args.replace is None:
		rep = input("Replacement string [default: AmgenrabbitmAb]: ").strip()
		replacement = rep or 'AmgenrabbitmAb'
	else:
		replacement = args.replace

	# Determine whether to perform renames now. If --apply was passed, perform
	# them; otherwise prompt the user (default No).
	if args.apply:
		dry_run = False
	else:
		ans = input('Perform renames now? (y/N): ').strip().lower()
		dry_run = False if ans in ('y', 'yes') else True

	try:
		rename_in_folder(folder, pattern=pattern, replacement=replacement, dry_run=dry_run)
	except Exception as e:
		print('Error:', e)
		sys.exit(1)



