class Molecule:
    def __init__(self, pos_x, pos_y, static):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.static = static

class Direction_movement:
    def __init__(self, u_prob, r_prob, d_prob, l_prob):
        self.u_prob = u_prob
        self.r_prob = r_prob
        self.d_prob = d_prob
        self.l_prob = l_prob

class Movement_chance:
    def __init__(self, direction,prob):
        self.dir = direction
        self.prob = prob