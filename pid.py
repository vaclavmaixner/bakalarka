import math
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

class PID:
    def __init__(self, kp, ki, kd, derivator, integrator):
        self.kp=kp
        self.ki=ki
        self.kd=kd
        self.derivator=derivator
        self.integrator=integrator

        self.set_point = 2.0
        self.error = 0.0


    def update(self, current_value):
        self.error = self.set_point - current_value

        self.p_value = self.kp * self.error
        self.d_value = self.kd * (self.error - self.derivator)
        self.derivator = self.error

        self.integrator += self.error

        self.i_value = self.integrator * self.ki

        pid = self.p_value + self.i_value + self.d_value

        return pid


    def set_point(self, set_point):
        self.set_point = set_point
        self.integrator = 0.0
        self.derivator = 0.0

new_pid = PID(0.0, 0.0, 0.0, 0.0, 0.0)

# z = 100.0
# for i in range(0,10):
#     output=new_pid.update(z)
#     z+=output
#     #plt.plot(i,z)

#     plt.scatter(i,z)
#     print(output)

# plt.show()
    


#set_point(2)


# for i in range(0,100):
#     x = math.cos(i)
#     new_position = new_pid.update(x)
#     distance = new_position - x

#     plt.ion() #dela interaktivni prostredi
#     plt.plot(i, x, marker='x', color='b', linestyle='-')
#     plt.show()
#     print(new_position)

# plt.show()
# plt.waitforbuttonpress()
