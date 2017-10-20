#!/usr/bin/python

import random
import math
from time import sleep
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from timeit import default_timer as timer
import subprocess

#strukturyhumans
import molecule as molecule_import
import pid as pid_import
#import graphics as graphics_import
#import plot as myplt

NO_MOLECULES = 20
HEIGHT = 21
WIDTH = 21
ITERATIONS = HEIGHT*WIDTH
WAIT_TIME = 0.0000000000000001


list_of_molecules = []
direction_count_list = []

#vyrobi nahodne molekuly a nacpe je do

def setup_molecules(no_molecules):
    for k in range(0, no_molecules):
        pos_x = random.randint(0, HEIGHT)
        pos_y = random.randint(0, WIDTH)
        molecule = molecule_import.Molecule(k, pos_x, pos_y)
        list_of_molecules.append(molecule)

def setup_test_molecules():
    molecule1 = molecule_import.Molecule(1, 0, 0)
    molecule2 = molecule_import.Molecule(2, 0, 0)
    molecule3 = molecule_import.Molecule(3, HEIGHT, WIDTH)
    list_of_molecules.append(molecule1)
    list_of_molecules.append(molecule2)
    list_of_molecules.append(molecule3)


#vybere index molekuly nahodne ze seznamu vsech molekul
def choose_molecule():
    k = random.randint(0, len(list_of_molecules))
    #print("Chosen molecule is of index:", k)
    return k

def countOverlap(difference_x, difference_y):
    #prechod od pocitani v poli do pocitani v kartezske soustave souradnic
    dy = difference_x
    dx = difference_y
    #print("we got in countoverlap",difference_x,"and ", dy)
    overlap = 0

    if (dx == 1 or dx == -1):
        if (dy == 2):
            overlap = 1
        elif (dy == 1):
            overlap = 2
        elif (dy == 0):
            overlap = 3
        elif (dy == -1):
            overlap = 2
        elif (dy == -2):
            overlap = 1
    elif (dx == 0):
        if (dy == 2):
            overlap = 2
        elif (dy == 1):
            overlap = 4
        elif (dy == 0):
            overlap = 6
        elif (dy == -1):
            overlap = 4
        elif (dy == -2):
            overlap = 2
    return overlap

def checkNeighbours(index_of_chosen_mol, x, y):
    sum_of_overlap = 0
    for index in range(0,len(list_of_molecules)):
        # difference_x=0
        # difference_y=0
        if (index != index_of_chosen_mol):
            potential_neighbour = list_of_molecules[index]
            difference_x = x - (potential_neighbour.pos_x % WIDTH)
            difference_y = y - (potential_neighbour.pos_y % HEIGHT)
            #print(difference_x, "mezera", difference_y)

            # diff_x=-1 znamena, ze sousedici molekula sedi
            # nad nami vysetrovanou molekulou, y=-1 je napravo

            sum_of_overlap += countOverlap(difference_x,difference_y)
    return sum_of_overlap

def virtualMove(index_of_chosen_mol, direction):
    chosen_mol = list_of_molecules[index_of_chosen_mol]
    moved_mol = molecule_import.Molecule(1, 0, 0)
    #print(moved_mol.pos_x, moved_mol.pos_y)

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
    #print("moved molecule", moved_mol.pos_x, " a y ", moved_mol.pos_y)

    chance = checkNeighbours(index_of_chosen_mol, moved_mol.pos_x, moved_mol.pos_y)
    move_chance = molecule_import.Movement_chance(direction,chance)
    #print("kontrola move_chance ma hodnotu ", move_chance.dir, " a sance je", move_chance.prob)
    return move_chance

def convertOverlapToProbability(overlapAfterMove, overlapInitial, direction):
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

    # print(overlapAfterMove, "and the initial overlap is ", overlapInitial)
    delta_energy = (overlapInitial - overlapAfterMove)*c
    #print("DELTA ENERGY IS", delta_energy)
    energy_act = ENERGY_EQ + delta_energy/2 + (delta_energy*delta_energy)/(16*ENERGY_EQ)
    probability = diffusivity_coefficient*GAMMA*math.exp((-1)*(energy_act)/(BOLTZMANN*TEMPERATURE))
    #print("probability IS: ", probability)

    return probability


def moveChance(index_of_chosen_mol):
    chance_sum = 0

    #purge the previous direction count list
    #for i in range(0,len(direction_count_list)):
    #    del direction_count_list[i]


    #print("the chosen molecule is of index", index_of_chosen_mol)
    for direction in range(1,5):
        move_chance = virtualMove(index_of_chosen_mol, direction)
        #print("And the chance isssssssss: ", checkNeighbours(index_of_chosen_mol,0,0))
        convertedProbability = convertOverlapToProbability(move_chance.prob, checkNeighbours(index_of_chosen_mol,list_of_molecules[index_of_chosen_mol].pos_x,list_of_molecules[index_of_chosen_mol].pos_y), direction)
        move_chance.prob = convertedProbability

        chance_sum += convertedProbability
        #print ("the chance to move for direction",move_chance.dir,"is ", move_chance.prob)

        current_dir_move = molecule_import.Direction_movement(index_of_chosen_mol, move_chance.dir, move_chance.prob)
        direction_count_list.append(current_dir_move)

    # print("the total probability is ", chance_sum)

def moveChanceForAll():
    for i in range(0,len(list_of_molecules)):
        moveChance(i)


def printDirectionCountList():
    for i in range(0,len(direction_count_list)):    #, end=""
        print("The index is " + str(direction_count_list[i].index) + ", the direction is ")
        print(str(direction_count_list[i].dir) + ", the chance is " + str(direction_count_list[i].prob))

def graphics(iterations):
    # fig = plt.figure()
    plt.clf()
    OFFSET = 1

    plt.ion()
    plt.axis([0 - OFFSET, WIDTH + OFFSET, 0 - OFFSET - 1, HEIGHT + OFFSET - 1])

    currentAxis = plt.gca()


    # ax = fig.add_subplot(111)

    for index in range(0,len(list_of_molecules)):
        x = list_of_molecules[index].pos_x
        y = list_of_molecules[index].pos_y

        #jmeno molekuly v jejim prostredku
        #currentAxis.annotate(str(index), xy=(4, 0), xytext=(x + 1.5, y-1.1))
        currentAxis.annotate(str(index), xy=(4, 0), xytext=(x+0.4, y-0.6))

        #pro cyklicke okrajove podminky zavedeme modulo souradnice pro kazdy element bunky ve stylu:
        # 1 2 3
        # 4 5 6
        x1 = x % WIDTH
        y1 = y % HEIGHT
        x2 = (x + 1) % WIDTH
        y2 = y % HEIGHT
        x3 = (x + 2) % WIDTH
        y3 = y % HEIGHT

        x4 = x % WIDTH
        y4 = (y-1) % HEIGHT
        x5 = (x + 1) % WIDTH
        y5 = (y-1) % HEIGHT
        x6 = (x + 2) % WIDTH
        y6 = (y-1) % HEIGHT

        currentAxis.add_patch(Rectangle((x1, y1), 1, -1, alpha=0.4, facecolor="red",edgecolor="black"))
        currentAxis.add_patch(Rectangle((x2, y2), 1, -1, alpha=0.4, facecolor="red",edgecolor="black"))
        currentAxis.add_patch(Rectangle((x3, y3), 1, -1, alpha=0.4, facecolor="red",edgecolor="black"))
        currentAxis.add_patch(Rectangle((x4, y4), 1, -1, alpha=0.4, facecolor="red",edgecolor="black"))
        currentAxis.add_patch(Rectangle((x5, y5), 1, -1, alpha=0.4, facecolor="red",edgecolor="black"))
        currentAxis.add_patch(Rectangle((x6, y6), 1, -1, alpha=0.4, facecolor="red",edgecolor="black"))

        plt.plot([0, 0], [-1, HEIGHT-1], color='k', linestyle='-', linewidth=0.5)
        plt.plot([0, WIDTH], [-1, -1], color='k', linestyle='-', linewidth=0.5)
        plt.plot([WIDTH, WIDTH], [-1, HEIGHT-1], color='k', linestyle='-', linewidth=0.5)
        plt.plot([0, WIDTH], [HEIGHT-1, HEIGHT-1], color='k', linestyle='-', linewidth=0.5)

    plt.show()
    # print(delta_time)
    plt.pause(WAIT_TIME)

def pseudoGraphics():
    surface = [["0" for i in range(HEIGHT+1)] for j in range(WIDTH+1)]

    for k in range(0,len(list_of_molecules)):
        surface[list_of_molecules[k].pos_x][list_of_molecules[k].pos_y] = "X"

    for i in range(0,HEIGHT+1):
        for j in range(0,WIDTH+1):  #, end=""
            print(surface[i][j]," ")
        print(" ")

def check_cyclic_boundary_conditions():
    for i in range(0,len(list_of_molecules)):
        #print("Before", list_of_molecules[i].pos_x, " and ", list_of_molecules[i].pos_y)
        list_of_molecules[i].pos_x = (list_of_molecules[i].pos_x) % WIDTH
        list_of_molecules[i].pos_y = (list_of_molecules[i].pos_y) % HEIGHT
        #print("After", list_of_molecules[i].pos_x, " and ", list_of_molecules[i].pos_y)

def add_deposition_chance():
    deposition = molecule_import.Direction_movement((len(direction_count_list) + 1), 0, 0.1)
    direction_count_list.append(deposition)

time_list = []

whole_time = 0
def print_time(delta_time):
    global whole_time
    whole_time += delta_time

    molecule_over_time = molecule_import.MoleculesOverTime(whole_time, len(list_of_molecules))
    time_list.append(molecule_over_time)

    # line = (whole_time, ' ', len(list_of_molecules))
    # with open('molecules_over_time_data4.txt', 'a') as f:
    #     f.write(str(whole_time))
    #     f.write(' ')
    #     f.write(str(len(list_of_molecules)))
    #     f.write('\n')


delta_time=0.0
def makeOneChange():
    sumOfAllMoves = 0
    indexFinder = 0
    indexOfchosenEvent = 0

    #check_cyclic_boundary_conditions()




    #secte vsechny cestnosti dohromady
    for i in range(0,len(direction_count_list)):
        sumOfAllMoves += direction_count_list[i].prob
    #print("THE SUM OF ALL MOVES IS ",sumOfAllMoves)

    a = random.uniform(0.0, 1.0)
    # promenlivy cas, nastaveny tak, ze pro vysoce pravdebopodobnou udalost se to stane rychle a naopak
    delta_time = ((-1)*math.log(a)/sumOfAllMoves)

    print_time(delta_time)

    #print("THE TIME IS ", delta_time, "exponent is ", math.log10(delta_time))
    #print(delta_time * pow(10,-(math.log10(delta_time))))
    #sleep(delta_time * pow(10,-(math.log10(delta_time))) / 10)

    #vybere dany ukazatel na cetnosti, ktery urci, jak udalost se uskutecni
    chosenCount = random.uniform(0.0, sumOfAllMoves)

    #prochazi vsechny udalosti s jejich cetnostmi a vybere tu, ktera sedi s ukazatelem
    for i in range(0,len(direction_count_list)):
        indexFinder += direction_count_list[i].prob
        #print("jaa",indexFinder, " ", chosenCount)
        if indexFinder >= chosenCount:
            indexOfchosenEvent = i
            break

    #actually move
    if direction_count_list[indexOfchosenEvent].dir == 1:
        list_of_molecules[direction_count_list[indexOfchosenEvent].index].pos_y += 1
    elif direction_count_list[indexOfchosenEvent].dir == 2:
        list_of_molecules[direction_count_list[indexOfchosenEvent].index].pos_x += 1
    elif direction_count_list[indexOfchosenEvent].dir == 3:
        list_of_molecules[direction_count_list[indexOfchosenEvent].index].pos_y -= 1
    elif direction_count_list[indexOfchosenEvent].dir == 4:
        list_of_molecules[direction_count_list[indexOfchosenEvent].index].pos_x -= 1
    elif direction_count_list[indexOfchosenEvent].dir == 0:
        setup_molecules(1)



    #print(direction_count_list[indexOfchosenEvent].index, "and the direction is ", direction_count_list[indexOfchosenEvent].dir)
    return indexOfchosenEvent

def print_molecule_list():
    for l in range(0,len(list_of_molecules)):
        print('The positions of molecules are: ',str(list_of_molecules[l].pos_x) + " " + str(list_of_molecules[l].pos_y))

def tip_neighbours(x,y):
    overlap=0.0
    for l in range(0,len(list_of_molecules)):
        potential_neighbour = list_of_molecules[l]
        difference_x = x - (potential_neighbour.pos_x % WIDTH)
        difference_y = y - (potential_neighbour.pos_y % HEIGHT)
            #print(difference_x, "mezera", difference_y)

            # diff_x=-1 znamena, ze sousedici molekula sedi
            # nad nami vysetrovanou molekulou, y=-1 je napravo

        sum_of_overlap = countOverlap(difference_x,difference_y)
        if sum_of_overlap!=0:
            overlap+=1
    return overlap


def call_PID(i,z,new_pid):
    pos_width=i%HEIGHT
    pos_height=(i-pos_width)/HEIGHT
    #sem kouka hrot mikroskopu, vystupem je vyska baliku molekul pod hrotem
    stack_height=tip_neighbours(pos_width,pos_height)
    if stack_height>1:
        #print(stack_height)
        pass

    update=new_pid.update(z-stack_height)
    #print(update)
    return(update)

def gnuplot():
    proc = subprocess.Popen(['gnuplot --persist', '-p'], 
                            shell=True,
                            stdin=subprocess.PIPE,
                            )
    # proc.stdin.write("set dgrid3d\n")
    # proc.stdin.write("splot 'out.txt' u 1:2:3 with pm3d\n")

    proc.stdin.write("set pm3d map\n")
    proc.stdin.write("set palette rgbformulae 22,13,-31\n")
    proc.stdin.write("splot 'out.txt' u 1:2:3\n")

def Main():
    start = timer()
    setup_molecules(NO_MOLECULES)

    new_z=0.0
    z=5.0
    new_pid = pid_import.PID(1.2, 0.1, 0.1, 0.0, 0.0)



    #choose_molecule()
    #moveChanceForAll()
    #graphics()
    #pseudoGraphics()

    #print(makeOneChange(), " is the chosen index of event to be done")

    #printDirectionCountList()

    #convertOverlapToProbability(1)

    print("TEEEST ", checkNeighbours(2, 7, 7))


    



    meanwhile_t=0.0
    #the main loop
    keep_running = True
    with open('out.txt', 'wt') as out:
        for i in range(0,ITERATIONS):
            add_deposition_chance()
            moveChanceForAll()
            #pseudoGraphics()
            #graphics(i)
            #graphics_import.graphics_main(list_of_molecules)
            #printDirectionCountList()
            #print_molecule_list()
            makeOneChange()

            dt=60/(512*20*512)
            meanwhile_t+=delta_time
            if meanwhile_t>=dt:
                meanwhile_t=0.0
                z+=call_PID(i,z,new_pid)

                pos_width=i%HEIGHT
                pos_height=(i-pos_width)/HEIGHT

                out_str=str(pos_width)+"    "+str(pos_height)+"    "+str(z)+"\n"
                #if i>=WIDTH:    
                out.write(out_str)
                if ((i+1)%WIDTH==0):
                    out.write("\n")


            #print("another run_________________________")
            del direction_count_list[:]
    plt.show(block=True)

    for i in range(0,len(time_list)):
        with open('molecules_over_time_data5.txt', 'a') as f:
            f.write(str(time_list[i].time))
            f.write(' ')
            f.write(str(time_list[i].no_molecules))
            f.write('\n')


    gnuplot()

    end = timer()
    print(end-start, "seconds")


#run the setup
Main()






#udalost se vybira tak, ze ze vsech molekul a ze vsech jejich pohybu se vybere jeden pohyb u jedny molekuly,
#test - udelat misto repulze pritazlivou silu
#pozdeji udelat prepocitavani pozic jenom lokalne, aby se usetril cas
#urcit hranice, udela cyklicky

#difuzivita 0,6, 300 Kelvin, 10e13 1/s prefaktor divny gamma

#plus depozice s konstanti psti

#vypocitat si tabulku prekryvu na pravdepodobnost
#nastavit si difuzivitu


#Soucasny problem:
#Jaktoze nemame vztazenou cetnost na casovy krok? v soucasnosti je relativni
#v radech e+85, casovy krok v radech e-86. Je to spravne? Proc se v textu
#zminuje relativni cetnost pridruzena k casu? Jaky vyznam ma cas, vystupuje nekde
#ve vypoctech?
#Proc pro zmenenou ENERGY_EQ vychazi hodne rozdilne vysledky?

#vzgrafovat si depozici proti casu, nemelo by se to lisit s vetsim prekryvem
#zavest difuzivitu na energiich, ne listu
#zamrazit jednu molekulu pro stm

#512 pixelu, v kazdym 200 hodnot
#1/4 vteriny na radek
#elementarni kosticka je dejme tomu 4 angstromy
#obrazek typicky 10 nanometru
#udelat filtry na zpetnou vazbu