#!/usr/bin/python
from __future__ import division
import random
import math
import time
from time import sleep
from timeit import default_timer as timer
import subprocess
import pprint
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as MT
import matplotlib.cm as CM
from matplotlib.patches import Rectangle

import molecule3 as molecule_import
import pid as pid_import
import graphics as graphics_import

HEIGHT = 15
WIDTH = 15
ITERATIONS=HEIGHT*WIDTH/10

STEP_HEIGHT = 5
KP = 0.7

DT=60.0/(512*512)

list_of_molecules={}

def setup_molecules(no_molecules):
	for k in range(0, no_molecules):
		pos_x = random.randint(0, WIDTH)
		pos_y = random.randint(0, HEIGHT)
		molecule = molecule_import.Molecule(pos_x, pos_y, False,0,0,0,0)
		list_of_molecules[k]=molecule

def setup_molecules_test(no_molecules):
	for k in range(0, no_molecules):
		pos_x = WIDTH//2
		pos_y = HEIGHT//2
		molecule = molecule_import.Molecule(pos_x, pos_y, False,0,0,0,0)
		list_of_molecules[k]=molecule

def print_molecule_list():
    print('\n')
    for l in range(0,len(list_of_molecules)):
        print('The positions of molecule', l, ' are: ',str(list_of_molecules[l].pos_x) + " " + str(list_of_molecules[l].pos_y))
    print('\n')

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

#takes in index of molecule chosen to check neigbouors for, takes in its x a y position
def check_neighbours(index_of_chosen_mol, x, y):
	sum_of_overlap = 0

	for index in range(0,len(list_of_molecules)):
		if (index != index_of_chosen_mol):
			potential_neighbour = list_of_molecules[index]
			difference_x = (x - (potential_neighbour.pos_x))
			difference_y = (y - (potential_neighbour.pos_y))

			if difference_x == WIDTH or difference_x == -WIDTH:
				difference_x == 1

			if difference_y == HEIGHT or difference_y == -HEIGHT:
				difference_y == 1

			if difference_x == WIDTH-1 or difference_x == -(WIDTH-1):
				difference_x == 2
			
			if difference_y == HEIGHT-1 or difference_y == -(HEIGHT-1):
				difference_y == 2

			if abs(difference_x)<=3 and abs(difference_y)<=1:
				sum_of_overlap += count_overlap(difference_x,difference_y)

			#print '\t ',index, ' dif x ', difference_x, ' dif y ', difference_y,
	return sum_of_overlap




def virtual_move(index_of_chosen_mol, direction):
	chosen_mol = list_of_molecules[index_of_chosen_mol]
	moved_mol = molecule_import.Molecule(0, 0, False,0,0,0,0)

	overlap_initial = check_neighbours(index_of_chosen_mol, chosen_mol.pos_x, chosen_mol.pos_y)
	#print 'initial'
	if direction==1:
		moved_mol.pos_x = chosen_mol.pos_x
		moved_mol.pos_y = chosen_mol.pos_y - 1
	elif direction==2:
		moved_mol.pos_x = chosen_mol.pos_x + 1
		moved_mol.pos_y = chosen_mol.pos_y
	elif direction==3:
		moved_mol.pos_x = chosen_mol.pos_x
		moved_mol.pos_y = chosen_mol.pos_y + 1
	elif direction==4:
		moved_mol.pos_x = chosen_mol.pos_x - 1
		moved_mol.pos_y = chosen_mol.pos_y

	overlap_after_move = check_neighbours(index_of_chosen_mol, moved_mol.pos_x, moved_mol.pos_y)
	#print 'after'
	move_overlap_count = molecule_import.Movement_chance(direction,overlap_after_move)

	#prob_after_move = convert_overlap_to_prob(overlap_after_move,overlap_initial, direction)

	return move_overlap_count


def convert_overlap_to_prob(overlapAfterMove, overlapInitial, direction):
	BOLTZMANN = 8.617e-5 #prepsat do eV
	TEMPERATURE = 300
	ENERGY_EQ = 0.8
	c = -0.1

	#print overlapAfterMove, overlapInitial, direction

	GAMMA = 1e+13

	#difuzivni koeficient pro pohyb povrchem nahoru/dolu a do stran
	vertical_coefficient = 1
	horizontal_coefficient = 1

	if (direction == 1) or (direction == 3):
		diffusivity_coefficient = vertical_coefficient
	elif (direction == 2) or (direction == 4):
		diffusivity_coefficient = horizontal_coefficient


	#tady je prusvih, ze ta energie muze vyjit zaporne a to je proste blbost,
	# bude to mit desnou sanci se stat, co nejak urcit konstantu c? nema rozmer?
	# ***************************
	delta_energy = (overlapInitial - overlapAfterMove)*c
	#print delta_energy, 'delta_energy		',

	energy_act = ENERGY_EQ + delta_energy/2 + (delta_energy*delta_energy)/(16*ENERGY_EQ)
	#print energy_act
	#print energy_act, ' energy act'
	probability = diffusivity_coefficient*GAMMA*math.exp((-1)*(energy_act)/(BOLTZMANN*TEMPERATURE))

	print probability, 'probability pro ', overlapAfterMove, "overlap after and ", overlapInitial, ' overlapInitial'
	return probability

def count_movement_chances(index_of_chosen_mol):
	chance_sum=0
	u_prob=0
	r_prob=0
	d_prob=0
	l_prob=0

	chosen_mol = list_of_molecules[index_of_chosen_mol]
	overlap_initial = check_neighbours(index_of_chosen_mol, chosen_mol.pos_x, chosen_mol.pos_y)
	#print 'initial2'
	for direction in range(1,5):
		overlap_after_move = virtual_move(index_of_chosen_mol, direction)
		print index_of_chosen_mol, 'overlap ', overlap_after_move.prob, ' initial ', overlap_initial, '\n'

		chance_after_move=convert_overlap_to_prob(overlap_after_move.prob, overlap_initial, direction)

		chance_sum+=chance_after_move
		if direction==1:
			list_of_molecules[index_of_chosen_mol].u_prob=chance_after_move
		elif direction==2:
			list_of_molecules[index_of_chosen_mol].r_prob=chance_after_move
		elif direction==3:
			list_of_molecules[index_of_chosen_mol].d_prob=chance_after_move
		elif direction==4:
			list_of_molecules[index_of_chosen_mol].l_prob=chance_after_move



def setup_directions():
	for i in range(0,len(list_of_molecules)):
		if not list_of_molecules[i].static:
			count_movement_chances(i)

def make_one_change():
	#print ".",
	sumOfAllMoves = 0
	indexFinder = 0
	indexOfchosenEvent = 0
	chosen_direction=0
	partial_sum = {}
	lil_sum=0.0

	for i in range(0,len(list_of_molecules)):
		
		#lil_sum=list_of_molecules[i].u_prob + list_of_molecules[i].r_prob + list_of_molecules[i].d_prob + list_of_molecules[i].l_prob
		lil_sum+=list_of_molecules[i].u_prob
		partial_sum[((i+1)*4)-3]=lil_sum
		#print lil_sum, 
		lil_sum+=list_of_molecules[i].r_prob
		partial_sum[((i+1)*4)-2]=lil_sum
		#print lil_sum, 
		lil_sum+=list_of_molecules[i].d_prob
		partial_sum[((i+1)*4)-1]=lil_sum
		#print lil_sum, 
		lil_sum+=list_of_molecules[i].l_prob
		partial_sum[((i+1)*4)]=lil_sum
		#print lil_sum
		print(i,'Probabilities:',list_of_molecules[i].u_prob, list_of_molecules[i].r_prob, list_of_molecules[i].d_prob, list_of_molecules[i].l_prob)


	sumOfAllMoves += lil_sum
	#print 'sum of all moves', sumOfAllMoves

		
	
	for i in range(1,len(partial_sum)+1):
		#print i, partial_sum[i], ' sumaa'
		pass

	a = random.uniform(0.0, 1.0)
	# promenlivy cas, nastaveny tak, ze pro vysoce pravdebopodobnou udalost se to stane rychle a naopak
	#delta_time = ((-1)*math.log(a)/sumOfAllMoves)

	chosenCount = random.uniform(0, sumOfAllMoves)
	print 'sum of all moves ', sumOfAllMoves, ' chosen count ', chosenCount
	#print('chosenCount')
	#print(sumOfAllMoves)
	#print(chosenCount)
	#print(' chosenCount\n')

	#prochazi vsechny udalosti s jejich cetnostmi a vybere tu, ktera sedi s ukazatelem
	

	for i in range (1,len(partial_sum)+1):
		if i==1:
			previous=0
		else:
			previous=partial_sum[i-1]
		if chosenCount<=partial_sum[i] and chosenCount>=previous:
			indexOfchosenEvent=i//4
			if i%4==0:
				indexOfchosenEvent-=1
			chosen_direction = i%4
			if chosen_direction==0:
				chosen_direction=4

			print indexOfchosenEvent, ' index a direction je', chosen_direction

	#actually move
	if chosen_direction == 1:
		list_of_molecules[indexOfchosenEvent].pos_y = (list_of_molecules[indexOfchosenEvent].pos_y - 1)%HEIGHT
	elif chosen_direction  == 2:
		 list_of_molecules[indexOfchosenEvent].pos_x = (list_of_molecules[indexOfchosenEvent].pos_x + 1)%WIDTH
	elif chosen_direction  == 3:
		 list_of_molecules[indexOfchosenEvent].pos_y = (list_of_molecules[indexOfchosenEvent].pos_y + 1)%HEIGHT
	elif chosen_direction  == 4:
		 list_of_molecules[indexOfchosenEvent].pos_x = (list_of_molecules[indexOfchosenEvent].pos_x - 1)%WIDTH

	print(indexOfchosenEvent, chosen_direction, 'pohyb molekuly')
	return indexOfchosenEvent


def loop():
	NO_MOLECULES = 10
	#setup_molecules(NO_MOLECULES)
	setup_molecules_test(NO_MOLECULES)
	print_molecule_list()

	#print_molecule_list()
	setup_directions()
	#print(list_of_molecules[0].u_prob)
	#print(list_of_molecules[0].r_prob)
	#print(list_of_molecules[0].d_prob)
	#print(list_of_molecules[0].l_prob)

	#print(check_neighbours(0, list_of_molecules[0].pos_x, list_of_molecules[0].pos_y)), 'puvodni'
	#print(virtual_move(0,1).prob, convert_overlap_to_prob(virtual_move(0,1).prob, check_neighbours(0, list_of_molecules[0].pos_x, list_of_molecules[0].pos_x), 1))
	#print(virtual_move(0,2).prob, convert_overlap_to_prob(virtual_move(0,2).prob, check_neighbours(0, list_of_molecules[0].pos_x, list_of_molecules[0].pos_x), 2))
	#print(virtual_move(0,3).prob, convert_overlap_to_prob(virtual_move(0,3).prob, check_neighbours(0, list_of_molecules[0].pos_x, list_of_molecules[0].pos_x), 3))
	#print(virtual_move(0,4).prob, convert_overlap_to_prob(virtual_move(0,4).prob, check_neighbours(0, list_of_molecules[0].pos_x, list_of_molecules[0].pos_x), 4))
	print('\n')

	for i in range(0,200):
		
		make_one_change()
		graphics_import.graphics(WIDTH, HEIGHT, list_of_molecules)
		setup_directions()
		print_molecule_list()
		print "AAAAAAAAAAAAAAAAAAAAAAAA"

loop()

#convert_overlap_to_prob(20,5,2)
#convert_overlap_to_prob(13,20,2) #tohle by melo byt hodne pravdepodobne, nejvice
#convert_overlap_to_prob(14,20,2)
#convert_overlap_to_prob(15,20,2)
#convert_overlap_to_prob(16,20,2)
#convert_overlap_to_prob(17,20,2) #maximum, takze to je dost jety
#convert_overlap_to_prob(18,20,2) #maximum, takze to je dost jety