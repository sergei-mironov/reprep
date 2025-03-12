#!/usr/bin/env python

import subprocess
import argparse
from os import walk, chdir
from os.path import isdir, isfile, join, basename

def search(entity, file_name) -> int|None:
  line_number = None
  with open(file_name, 'r') as file:
    for i, line in enumerate(file, start=1):
      if f"class {entity}" in line.strip():
        line_number = i
        break
  return line_number


def scan(url_template, entities, file_path):
  result = []
  pyfiles = []
  # Check if the file_path is a directory
  if isdir(file_path):
    for root, dirs, files in walk(file_path):
      # Filter-out directories with dot-names (e.g., .git, .venv)
      dirs[:] = [d for d in dirs if not d.startswith('.')]
      # Loop through each file in the directory
      for file in files:
        # Consider only .py files
        if file.endswith('.py'):
          full_file_path = join(root, file)
          pyfiles.append(full_file_path)
  else:
    pyfiles.append(file_path)

  for e in entities:
    if isfile(e):
      result.append(url_template
                    .replace('%E',basename(e))
                    .replace('%F',e)
                    .replace('%L','1'))
    else:
      for f in pyfiles:
        ln = search(e, f)
        if ln is not None:
          result.append(url_template
                        .replace('%E',e)
                        .replace('%F',f)
                        .replace('%L',str(ln)))
  return result

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description='Generate links for class references in a Python file.')
  parser.add_argument('entities', metavar='E', type=str, nargs='+',
                      help='entities to generate links for')
  parser.add_argument('--path', metavar='F', type=str, default='.',
                      help='Path to the Python file to search for classes')
  parser.add_argument('-C', metavar='DIR', type=str, default=None,
                      help='set the working directory before execution')
  parser.add_argument('--url', type=str, default="\\href{%R/blob/%B/%F#L%L}{%E}",
                      help=('URL template for generating class links; %%R - remote; %%B - git '
                      'revision; %%F - file name; %%L - line number; %%E - entity name'))
  parser.add_argument('-S', '--separator', type=str, default=' | ',
                      help='Separator used to join the URL links in the output')

  args = parser.parse_args()

  if args.C is None:
    args.C = '.'
  else:
    chdir(args.C)
    if args.path.startswith(args.C+'/'):
      args.path = args.path[len(args.C+'/'):]

  # Calculate the git_remote by calling the system's git to get the remote URL
  git_remote = subprocess.check_output(['git', 'config', '--get', 'remote.origin.url'],
                                       text=True).strip()

  # Calculate the git_rev by calling the system's git to get the current HEAD revision
  git_rev = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()

  # Replace placeholders in the URL template
  url_template = args.url.replace('%C',args.C).replace('%R', git_remote).replace('%B', git_rev)
  result = scan(url_template, args.entities, args.path)
  print(args.separator.join(result))

