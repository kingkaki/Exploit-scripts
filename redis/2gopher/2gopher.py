# -*- coding: utf-8 -*-
# @Author: kingkk
# @Date:   2018-08-18 10:33:18
# @Last Modified by:   kingkk
# @Last Modified time: 2018-08-18 12:51:18
import sys

exp = ''

with open(sys.argv[1],'r') as f:
	for line in f.readlines():
		if line[0] in '><+':
			continue
		elif line[-3:-1] == r'\r':
			if len(line) == 3:
				exp = exp + '%0a%0d%0a'
			else:
				line = line.replace(r'\r', '%0d%0a')
				line = line.replace('\n', '')
				exp = exp + line
		elif line == '\x0a':
			exp = exp + '%0a'
		else:
			line = line.replace('\n', '')
			exp = exp + line
print(exp)