from flask import *
import pandas as pd
import numpy as np
from fileinput import filename
from werkzeug.utils import secure_filename
import os

UPLOAD_FOLDER = '/flask_hydrogene/uploads'
ALLOWED_EXTENSIONS = {'out'}


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['SECRET_KEY'] = 'string'

@app.route("/")
def home():
	return render_template('index.html')

@app.route("/uploaded", methods = ['GET', 'POST'])	
def uploaded() -> "html":
	if request.method == 'POST':
		f=request.files['file']
		f.save(f.filename)
		with open(f.filename) as out_file:
			lines = out_file.readlines()
		n_atom_1 = request.form['n_atom_1']
		n_atom_1 = int(n_atom_1)
		n_atoms_start = 'xyz'
		n_atoms_end = 'end of input'
		xyz_list1_h, xyz_list2_h, xyz_list1_h_index, xyz_list2_h_index = [], [], [], []
		xyz_list_cycle1_h, xyz_list_cycle2_h, dist_init_list, dist_cycle_list = [], [], [], []
		for l_index, line in enumerate(lines):
			line = line.lower().strip()
			if n_atoms_start in line:
				l_index_start = l_index
				l_index_start = l_index_start
			elif n_atoms_end in line:
				l_index_end = l_index
				l_index_end = l_index_end
		num_atoms = l_index_end - l_index_start -1
		print(f"number of atoms: {num_atoms}, molecule 1 - {n_atom_1} atoms, molecule 2 - {num_atoms - n_atom_1} atoms")
		xyz_list = lines[l_index_start+1:l_index_end]
		xyz_list1 = xyz_list[:n_atom_1]										
		xyz_list2 = xyz_list[n_atom_1:]
		for item in xyz_list1:
			item = item.replace('|', ' ').replace('>', ' ')
			if item.find('H') != -1:
				xyz_list1_h.append(item.strip().split())
		for item in xyz_list2:
			item = item.replace('|', ' ').replace('>', ' ')
			if item.find('H') != -1:
				xyz_list2_h.append(item.strip().split())
		for i in range(len(xyz_list1_h)):
			item = xyz_list1_h[i][0]
			xyz_list1_h_index.append(item)
		for k in range(len(xyz_list2_h)):
			item = xyz_list2_h[k][0]
			xyz_list2_h_index.append(item)
		xyz_h1_df = pd.DataFrame(xyz_list1_h, columns = ['atom number', 'atom', 'x', 'y', 'z'])
		xyz_h2_df = pd.DataFrame(xyz_list2_h, columns = ['atom number', 'atom', 'x', 'y', 'z'])
		xyz_h1_df.loc[:, ['x','y','z']] = xyz_h1_df.loc[:, ['x','y','z']].astype(float)
		xyz_h2_df.loc[:, ['x','y','z']] = xyz_h2_df.loc[:, ['x','y','z']].astype(float)
		for i in range(len(xyz_list1_h)):
			a = np.array([xyz_h1_df.iloc[i,2],xyz_h1_df.iloc[i,3],xyz_h1_df.iloc[i,4]])
			for k in range(len(xyz_list2_h)):
				b = np.array([xyz_h2_df.iloc[k,2],xyz_h2_df.iloc[k,3],xyz_h2_df.iloc[k,4]])
				dist = np.linalg.norm(a-b)
				if dist < 3.0:
					dist_init = f"Distance between: {xyz_list1_h[i][0]} and: {xyz_list2_h[k][0]} is: {dist}"
					dist_init_list.append(dist_init) 
					print(f"Initial Geometry: Distance between {xyz_list1_h[i][0]} and {xyz_list2_h[k][0]} atom: {dist}nm")
		opt_lines = 0
		for index, line in enumerate(lines):
			line = line.lower().strip()
			if 'geometry optimization cycle' in line:
				opt_lines += 1
		print(f" Number of optimization cycles: {opt_lines}")
		cycle_n = request.form['cycle_n']
		cycle_n = int(cycle_n)
		if cycle_n < 10:
			cycle_string = f"cycle   {cycle_n}"
		elif 10<= cycle_n <100:
			cycle_string = f"cycle  {cycle_n}"
		else:
			cycle_string = f"cycle {cycle_n}"
		res_dist = request.form['res_dist']
		res_dist = float(res_dist)											
		for index, line in enumerate(lines):
			line = line.lower().strip()
			if cycle_string in line:
				s_index = index +6
				e_index = s_index + num_atoms
				xyz_list_cycle = lines[s_index:e_index]
				xyz_list_cycle1 = xyz_list_cycle[:n_atom_1]
				xyz_list_cycle2 = xyz_list_cycle[n_atom_1:]
		for item in xyz_list_cycle1:
			item = item.replace('|', ' ').replace('>', ' ')
			if item.find('H') != -1:
				xyz_list_cycle1_h.append(item.strip().split())
		for item in xyz_list_cycle2:
			item = item.replace('|', ' ').replace('>', ' ')
			if item.find('H') != -1:
				xyz_list_cycle2_h.append(item.strip().split())
		xyz_cycle1_df = pd.DataFrame(xyz_list_cycle1_h, columns = ['atom', 'x', 'y', 'z'])
		xyz_cycle2_df = pd.DataFrame(xyz_list_cycle2_h, columns = ['atom', 'x', 'y', 'z'])
		xyz_cycle1_df['atom number'] = xyz_list1_h_index
		xyz_cycle2_df['atom number'] = xyz_list2_h_index
		xyz_cycle1_df = xyz_cycle1_df[['atom number', 'atom', 'x', 'y', 'z']]
		xyz_cycle2_df = xyz_cycle2_df[['atom number', 'atom', 'x', 'y', 'z']]
		xyz_cycle1_df.loc[:, ['x','y','z']] = xyz_cycle1_df.loc[:, ['x','y','z']].astype(float)
		xyz_cycle2_df.loc[:, ['x','y','z']] = xyz_cycle2_df.loc[:, ['x','y','z']].astype(float)
		for j in range(len(xyz_list_cycle1_h)):
			c = np.array([xyz_cycle1_df.iloc[j,2],xyz_cycle1_df.iloc[j,3],xyz_cycle1_df.iloc[j,4]])
			for l in range(len(xyz_list_cycle2_h)):
				d = np.array([xyz_cycle2_df.iloc[l,2],xyz_cycle2_df.iloc[l,3],xyz_cycle2_df.iloc[l,4]])
				dist = np.linalg.norm(c-d)
				if dist < res_dist:
					print(f"Cycle {cycle_n}: Distance between {xyz_cycle1_df.iloc[j,0]} and {xyz_cycle2_df.iloc[l,0]} atom: {dist}nm")
					dist_cycle = f"Cycle {cycle_n}: Distance between {xyz_cycle1_df.iloc[j,0]} and {xyz_cycle2_df.iloc[l,0]} atom: {dist}nm"
					dist_cycle_list.append(dist_cycle)
		return render_template('uploaded.html', fname = f.filename, num_atoms=num_atoms, opt_lines=opt_lines, dist_init_list=dist_init_list, dist_cycle_list=dist_cycle_list)


if __name__ == "__main__":
	app.run(debug = True)