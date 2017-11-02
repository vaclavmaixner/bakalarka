#!/usr/bin/python

import random
import math
from time import sleep
from timeit import default_timer as timer
import subprocess
import pprint
import profile
import re
import cProfile

NO_MOLECULES = 100
HEIGHT = 10
WIDTH = 20*10
ITERATIONS=HEIGHT*WIDTH

DT=60.0/(512*20*512)

list_of_molecules={}
directions_list={}

import molecule2 as molecule_import
import pid as pid_import

def setup_molecules(no_molecules):
	for k in range(0, no_molecules):
		pos_x = random.randint(0, WIDTH)
		pos_y = random.randint(0, HEIGHT)
		molecule = molecule_import.Molecule(pos_x, pos_y, False)
		list_of_molecules[k]=molecule

def append_static_molecule(x,y):
	k=len(directions_list)+1
	molecule=molecule_import.Molecule(pos_x, pos_y, True)
	directions_list[k]=molecule


def count_overlap(dx, dy):
	overlap = 0

	if (dy == 1 or dy == -1):
		if (dx == 2):
			overlap = 1
		elif (dx == 1):
			overlap = 2
		elif (dx == 0):
			overlap = 3
		elif (dx == -1):
			overlap = 2
		elif (dx == -2):
			overlap = 1
	elif (dy == 0):
		if (dx == 2):
			overlap = 2
		elif (dx == 1):
			overlap = 4
		elif (dx == 0):
			overlap = 6
		elif (dx == -1):
			overlap = 4
		elif (dx == -2):
			overlap = 2
	return overlap

def check_influential_neighbours(index_of_chosen_mol):
	list_of_inf_neighbours=set()

	x=list_of_molecules[index_of_chosen_mol].pos_x
	y=list_of_molecules[index_of_chosen_mol].pos_y

	for index in range(0,len(list_of_molecules)):
		if (index != index_of_chosen_mol):
			potential_neighbour = list_of_molecules[index]
			difference_x = x - (potential_neighbour.pos_x % WIDTH)
			difference_y = y - (potential_neighbour.pos_y % HEIGHT)

			if abs(difference_x)<=3 and abs(difference_y)<=1:
				list_of_inf_neighbours.add(index)

	return list_of_inf_neighbours

def check_neighbours(index_of_chosen_mol, x, y):
	sum_of_overlap = 0

	for index in range(0,len(list_of_molecules)):
		if (index != index_of_chosen_mol):
			potential_neighbour = list_of_molecules[index]
			difference_x = x - (potential_neighbour.pos_x % WIDTH)
			difference_y = y - (potential_neighbour.pos_y % HEIGHT)

			if abs(difference_x)<=3 and abs(difference_y)<=1:
				sum_of_overlap += count_overlap(difference_x,difference_y)
	return sum_of_overlap

def virtual_move(index_of_chosen_mol, direction):
	chosen_mol = list_of_molecules[index_of_chosen_mol]
	moved_mol = molecule_import.Molecule(0, 0, False)

	if direction==1:
		moved_mol.pos_x = chosen_mol.pos_x
		moved_mol.pos_y = chosen_mol.pos_y + 1
	elif direction==2:
		moved_mol.pos_x = chosen_mol.pos_x + 1
		moved_mol.pos_y = chosen_mol.pos_y
	elif direction==3:
		moved_mol.pos_x = chosen_mol.pos_x
		moved_mol.pos_y = chosen_mol.pos_y - 1
	elif direction==4:
		moved_mol.pos_x = chosen_mol.pos_x - 1
		moved_mol.pos_y = chosen_mol.pos_y

	overlap_after_move = check_neighbours(index_of_chosen_mol, moved_mol.pos_x, moved_mol.pos_y)
	move_overlap_count = molecule_import.Movement_chance(direction,overlap_after_move)

	return move_overlap_count


def convert_overlap_to_prob(overlapAfterMove, overlapInitial, direction):
	BOLTZMANN = 8.617e-5 #prepsat do eV
	TEMPERATURE = 300
	ENERGY_EQ = 0.8
	c = -1
	GAMMA = 1e+13

	#difuzivni koeficient pro pohyb povrchem nahoru/dolu a do stran
	vertical_coefficient = 1
	horizontal_coefficient = 2

	if (direction == 1) or (direction == 3):
		diffusivity_coefficient = vertical_coefficient
	elif (direction == 2) or (direction == 4):
		diffusivity_coefficient = horizontal_coefficient

	delta_energy = (overlapInitial - overlapAfterMove)*c
	energy_act = ENERGY_EQ + delta_energy/2 + (delta_energy*delta_energy)/(16*ENERGY_EQ)
	probability = diffusivity_coefficient*GAMMA*math.exp((-1)*(energy_act)/(BOLTZMANN*TEMPERATURE))

	return probability

#pridat boolean updat=True, pak je append, jinak menit jenom dane hondoty
def count_movement_chances(index_of_chosen_mol):
	chance_sum=0
	u_prob=0
	r_prob=0
	d_prob=0
	l_prob=0

	overlap_initial=check_neighbours

	for direction in range(1,5):
		overlap_after_move = virtual_move(index_of_chosen_mol, direction)

		x=list_of_molecules[index_of_chosen_mol].pos_x
		y=list_of_molecules[index_of_chosen_mol].pos_y

		overlap_initial=check_neighbours(index_of_chosen_mol,x,y)
		chance_after_move=convert_overlap_to_prob(overlap_after_move.prob, overlap_initial, direction)

		chance_sum+=chance_after_move
		if direction==1:
			u_prob=chance_after_move
		elif direction==2:
			r_prob=chance_after_move
		elif direction==3:
			d_prob=chance_after_move
		elif direction==4:
			l_prob=chance_after_move

	current_dir_move=molecule_import.Direction_movement(u_prob, r_prob, d_prob, l_prob)
	
	directions_list[index_of_chosen_mol]=current_dir_move

def setup_directions_list():
	for i in range(0,len(list_of_molecules)):
		if not list_of_molecules[i].static:
			count_movement_chances(i)



def make_one_change():
	sumOfAllMoves = 0
	indexFinder = 0
	indexOfchosenEvent = 0
	chosen_direction=0
	global time_change

	length=len(directions_list)

	for i in range(0,length):
		lil_sum=directions_list[i].u_prob + directions_list[i].r_prob + directions_list[i].d_prob + directions_list[i].l_prob
		sumOfAllMoves += lil_sum

	a = random.uniform(0.0, 1.0)
	# promenlivy cas, nastaveny tak, ze pro vysoce pravdebopodobnou udalost se to stane rychle a naopak
	delta_time = ((-1)*math.log(a)/sumOfAllMoves)

	print(delta_time)
	time_change=delta_time

	chosenCount = random.uniform(0.0, sumOfAllMoves)


	#prochazi vsechny udalosti s jejich cetnostmi a vybere tu, ktera sedi s ukazatelem
	for i in range(0,length):
		for j in range(1,5):
			if j==1:
				indexFinder+=directions_list[i].u_prob
				if indexFinder >= chosenCount:
					indexOfchosenEvent=i
					chosen_direction=j
					break
			elif j==2:
				indexFinder+=directions_list[i].r_prob
				if indexFinder >= chosenCount:
					indexOfchosenEvent=i
					chosen_direction=j
					break
			elif j==3:
				indexFinder+=directions_list[i].d_prob
				if indexFinder >= chosenCount:
					indexOfchosenEvent=i
					chosen_direction=j
					break
			elif j==4:
				indexFinder+=directions_list[i].l_prob
				if indexFinder >= chosenCount:
					indexOfchosenEvent=i
					chosen_direction=j
					break

	#actually move
	if chosen_direction == 1:
		list_of_molecules[indexOfchosenEvent].pos_y += 1
	elif chosen_direction  == 2:
		 list_of_molecules[indexOfchosenEvent].pos_x += 1
	elif chosen_direction  == 3:
		 list_of_molecules[indexOfchosenEvent].pos_y -= 1
	elif chosen_direction  == 4:
		 list_of_molecules[indexOfchosenEvent].pos_x -= 1

	#print(direction_count_list[indexOfchosenEvent].index, "and the direction is ", direction_count_list[indexOfchosenEvent].dir)
	return indexOfchosenEvent

time_change=0

def tip_neighbours(x,y):
	overlap=0.0
	for l in range(0,len(list_of_molecules)):
		potential_neighbour = list_of_molecules[l]
		difference_x = x - (potential_neighbour.pos_x % WIDTH)
		difference_y = y - (potential_neighbour.pos_y % HEIGHT)
		sum_of_overlap = count_overlap(difference_x,difference_y)

		if sum_of_overlap!=0:
			overlap+=1
			
	return overlap

def call_PID(pos_width,pos_height,z,new_pid):
	#sem kouka hrot mikroskopu, vystupem je vyska baliku molekul pod hrotem
	stack_height=tip_neighbours(pos_width,pos_height)

	update=new_pid.update(z-stack_height)
	return(update)

def gnuplot():
	proc = subprocess.Popen(['gnuplot --persist', '-p'], 
							shell=True,
							stdin=subprocess.PIPE,
							)
	proc.stdin.write("set pm3d map\n")
	proc.stdin.write("set palette rgbformulae 22,13,-31\n")
	proc.stdin.write("set grid lt 1 lc 'black'\n")
	proc.stdin.write("set xtics 20\n")
	proc.stdin.write("set ytics 1\n")
	proc.stdin.write("splot 'out.txt' u 1:2:3\n")


def Main():
	

	start=timer()
	setup_molecules(NO_MOLECULES)
	setup_directions_list()
	
	z=5.0
	new_pid = pid_import.PID(1.4, 0.1, 0.1, 0.0, 0.0)

	print("start")

  
	with open('output2.txt', 'wt') as out:
		for i in range(0,ITERATIONS):
			meanwhile_t=0.0

			pos_width=i%WIDTH
			pos_height=(i-pos_width)/WIDTH

			while meanwhile_t<=DT:
				index_of_moved_mol=make_one_change()
				neighbours_to_check=check_influential_neighbours(index_of_moved_mol)
				for i in neighbours_to_check:
					count_movement_chances(i)

				print(time_change)

				meanwhile_t+=(time_change*10e22)

				z+=call_PID(pos_width,pos_height,z,new_pid)
				
			out_str=str(pos_width)+"    "+str(pos_height)+"    "+str(z)+"\n"
			out.write(out_str)
			if ((i+1)%WIDTH==0):
					out.write("\n")

	gnuplot()
	end=timer()
	print(end-start, "seconds")

Main()