#!/usr/bin/python
import os
import sys
def list_files(startpath,file_flag):
    global s
    for root, dirs, files in os.walk(startpath):
	level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        s += '{}{}/\n'.format(indent, os.path.basename(root))
        subindent = ' ' * 4 * (level + 1)
	if file_flag == 1 and os.path.basename(root) not in not_show_file:
		for f in files:
        		s += '{}{}\n'.format(subindent, f)
flag = raw_input('Display all files? (y/n)')
if flag == 'yes' or flag == 'y':
	file_flag = 1
else: file_flag = 0
s = ''
not_show_file = ['scratch','README']
only_show = []
list_files(os.getcwd(),file_flag)
f_out = open('dir_map.dat','w')
f_out.write(s)
