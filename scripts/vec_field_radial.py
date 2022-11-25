#!/usr/bin/env python
import rospy
from geometry_msgs.msg import Twist, Polygon, Pose
from turtlesim.msg import Pose as PoseTurtle
from nav_msgs.msg import Odometry
from std_msgs.msg import Bool
from sensor_msgs.msg import LaserScan
from math import sqrt, atan2, exp, atan, cos, sin, acos, pi, asin, atan2
from tf.transformations import euler_from_quaternion, quaternion_from_euler
from time import sleep
from visualization_msgs.msg import Marker, MarkerArray
import tf
from tf2_msgs.msg import TFMessage
# from scipy.spatial.transform import Rotation
import numpy as np
import copy
import sys
import copy
from time import gmtime, strftime
import rospkg


"""
Universidade Federal de Minas Gerais (UFMG) - 2019
Laboraorio CORO
Instituto Tecnologico Vale (ITV)
Contact:
Adriano M. C. Rezende, <adrianomcr18@gmail.com>
"""


# Obtain the path to be followed by the robot
def callback_path(data):

    global path, W
    global has_path


    path = [[],[]]
    # Iterate over the message to get all of the points of the curve
    for k in range(len(data.points)):
        p = data.points[k]
        if (invert_direction):
            path[0] = [p.x] + path[0]
            path[1] = [p.y] + path[1]
        else:
            path[0].append(p.x)
            path[1].append(p.y)


    Npts = len(path[0])


    D = 1.0
    F = [[],[]]
    points0 = copy.deepcopy(path)
    A0 = [[0] for i in range(Npts)]
    for k in range(1,Npts-1,1):
        T = [path[0][k+1]-path[0][k-1],path[1][k+1]-path[1][k-1]];
        norma = sqrt(T[0]**2+T[1]**2)
        T = [T[0]/norma, T[1]/norma]
        # N = np.matrix([[0, -1], [1, 0]])*T;
        N = [-T[1], T[0]]
        # N = np.array(N)

        path[0].append(path[0][k] + D*N[0]);
        path[1].append(path[1][k] + D*N[1]);
        A0.append([D]);
        path[0].append(path[0][k] - D*N[0]);
        path[1].append(path[1][k] - D*N[1]);
        A0.append([-D]);

    Npts2 = Npts + 2*(Npts-2);
    F0 = [[0 for j in range(Npts2)] for i in range(Npts2)]
    for i in range(Npts2):
        for j in range(Npts2):
            F0[j][i] = f_func([path[0][i], path[1][i]], [path[0][j], path[1][j]])


    A = np.matrix(A0)
    F = np.matrix(F0)
    W = np.linalg.inv(F)*(A + 1)


    #Save points to file
    rospack = rospkg.RosPack()
    file_path = rospack.get_path("vector_field_control") + "/log/coputed_alpha.txt"
    FILE = open(file_path, "w")
    #path_t = np.array(path).transpose().tolist()
    for k in path[0]:
        FILE.write(str(k) + "\t")
    FILE.write("\n")
    for k in path[1]:
        FILE.write(str(k) + "\t")
    FILE.write("\n")
    for k in W.tolist():
        FILE.write(str(k[0]) + "\t")
    #FILE.write("\n")

    FILE.close()

    has_path = True



#Radial basis function
def f_func(qa, qb):

    f = sqrt((qa[0]-qb[0])**2 + (qa[1]-qb[1])**2)

    return f


#Alpha function
def alpha_func(x,y,points,W):

    Npts = len(points[0])
    S = 0.0
    for k in range(Npts):
        S = S + W[k]*f_func([x, y],[points[0][k], points[1][k]])
    a = 1.0-S

    return a



# Callback to get the pose of the robot from tf
def callback_tf(data):
	global pos, rpy

	for T in data.transforms:
		# Chose the transform of the robot
		if (T.child_frame_id == "base_link"):

			# Get the orientation
			x_q = T.transform.rotation.x
			y_q = T.transform.rotation.y
			z_q = T.transform.rotation.z
			w_q = T.transform.rotation.w
			euler = euler_from_quaternion([x_q, y_q, z_q, w_q])
			theta_n = euler[2]

			# Get the position
			pos[0] = T.transform.translation.x
			pos[1] = T.transform.translation.y
			pos[2] = T.transform.translation.z
			rpy = euler

	return
# ----------  ----------  ----------  ----------  ----------





# Callback to get the pose of the robot
def callback_pose(data):
    global pos, rpy

    # Get the position
    pos[0] = data.position.x
    pos[1] = data.position.y
    pos[2] = data.position.z

    # Get the orientation
    x_q = data.orientation.x
    y_q = data.orientation.y
    z_q = data.orientation.z
    w_q = data.orientation.w
    euler = euler_from_quaternion([x_q, y_q, z_q, w_q])
    rpy = euler

    return
# ----------  ----------  ----------  ----------  ----------





# Callback to get the pose from odometry data
def callback_odometry(data):
    global pos, rpy

    # Get the position
    pos[0] = data.pose.pose.position.x
    pos[1] = data.pose.pose.position.y
    pos[2] = data.pose.pose.position.z

    # Get the orientation
    x_q = data.pose.pose.orientation.x
    y_q = data.pose.pose.orientation.y
    z_q = data.pose.pose.orientation.z
    w_q = data.pose.pose.orientation.w
    euler = euler_from_quaternion([x_q, y_q, z_q, w_q])
    rpy = euler

    return
# ----------  ----------  ----------  ----------  ----------




# Callback to get the pose from turtlesim data
def callback_turtlesim_pose(data):
    global pos, rpy

    # Get the position
    pos[0] = data.x
    pos[1] = data.y

    # Get the orientation
    rpy = [0, 0, data.theta]

    return
# ----------  ----------  ----------  ----------  ----------




# Get the alpha function and its gradient
def get_alpha_and_grad(x,y):

    global path, W

    #Delta for numerical derivative
    delta = 0.01

    #Compute alpha
    alpha = alpha_func(x,y,path,W)
    # print "alpha = ", alpha

    #Compute the gradient of alpha
    alpha_dx = alpha_func(x+delta,y,path,W)
    alpha_dy = alpha_func(x,y+delta,path,W)
    grad_x = (alpha_dx-alpha)/delta
    grad_y = (alpha_dy-alpha)/delta

    return alpha, grad_x, grad_y





# Compute the vectro field that will guide the robot
def vec_field_alpha(pos):

    global invert_direction
    global vd
    global d

    # Get the robot's x and y position
    x = pos[0]
    y = pos[1]

    # Get the position of the control point "p_d"
    psi = rpy[2]
    x = x + d*cos(psi)
    y = y + d*sin(psi)

    # Compute the alpha function and its gradients
    [alpha, grad_x, grad_y] = get_alpha_and_grad(x,y)

    print("[x, y] = %.2f, %.2f" % (x, y))
    print("[alpha, grad_x, grad_y] = %.2f, %.2f, %.2f" % (alpha, grad_x, grad_y))
    print ("")    

    # Normalize the gradient
    norm_grad = sqrt(grad_x**2+grad_y**2)
    grad_x_bar = grad_x/norm_grad
    grad_y_bar = grad_y/norm_grad

    # Copute the convergent and circulation weights
    k_G = -(2/pi)*atan(kf*alpha)
    k_H = sqrt(1-k_G**2)


    # Compute the field's componnents
    Vx = vd*( k_G*grad_x_bar - k_H*grad_y_bar )
    Vy = vd*( k_G*grad_y_bar + k_H*grad_x_bar )

    return (Vx, Vy)
# ----------  ----------  ----------  ----------  ----------



# Function feedback linearization
def feedback_linearization(Ux, Uy):
	global d

	# Get the yaw angle
	psi = rpy[2]

	# Compute foward velocity and angular velocity
	VX = cos(psi) * Ux + sin(psi) * Uy
	WZ = (-sin(psi) / d) * Ux + (cos(psi) / d) * Uy

	return (VX, WZ)
# ----------  ----------  ----------  ----------  ----------




# Function to send a markers, representing the value of the field
def send_marker_to_rviz(pub_rviz, Vx, Vy):

	mark_ref = Marker()

	mark_ref.header.frame_id = "world"
	mark_ref.header.stamp = rospy.Time.now()
	mark_ref.id = 0
	mark_ref.type = mark_ref.ARROW
	mark_ref.action = mark_ref.ADD
	# Size of the marker
	mark_ref.scale.x = 1.5 * (Vy ** 2 + Vx ** 2) ** (0.5)
	mark_ref.scale.y = 0.08
	mark_ref.scale.z = 0.08
	# Collor and transparency
	mark_ref.color.a = 1.0
	mark_ref.color.r = 0.0
	mark_ref.color.g = 0.0
	mark_ref.color.b = 0.0
	# Position of the marker
	mark_ref.pose.position.x = pos[0]
	mark_ref.pose.position.y = pos[1]
	mark_ref.pose.position.z = pos[2]
	#Orientation of the marker
	quaternio = quaternion_from_euler(0, 0, atan2(Vy, Vx))
	mark_ref.pose.orientation.x = quaternio[0]
	mark_ref.pose.orientation.y = quaternio[1]
	mark_ref.pose.orientation.z = quaternio[2]
	mark_ref.pose.orientation.w = quaternio[3]

	# Publish marker
	pub_rviz_ref.publish(mark_ref)

	return

# ----------  ----------  ----------  ----------  ----------


def load_parameters():

    global vd, kf, d
    global a, b, cx, cy, gamma, invert_direction
    global invert_motion_flag
    global pose_topic_name, pose_topic_type, cmd_vel_topic_name


    # Obtain the parameters
    try:
        vd = float(rospy.get_param("/vector_field/vd"))
        kf = float(rospy.get_param("/vector_field/K_F"))
        d = float(rospy.get_param("/vector_field/d_feedback"))

        invert_direction = bool(rospy.get_param("/vector_field/invert_direction"))
        invert_motion_flag = rospy.get_param("/vector_field/invert_motion_flag")
        pose_topic_name = rospy.get_param("/vector_field/pose_topic_name")
        pose_topic_type = rospy.get_param("/vector_field/pose_topic_type")
        cmd_vel_topic_name = rospy.get_param("/vector_field/cmd_vel_topic_name")
        path_topic_name = rospy.get_param("/vector_field/path_topic_name")


        print("\n\33[92mParameters loaded:\33[0m")
        print("\33[94mvr: " + str(vd) +"\33[0m")
        print("\33[94mkf: " + str(kf) +"\33[0m")
        print("\33[94md: " + str(d) +"\33[0m")
        print("\33[94minvert_direction: " + str(invert_direction) +"\33[0m")
        print("\33[94minvert_motion_flag: " +  str(invert_motion_flag) +"\33[0m")
        print("\33[94mpose_topic_name: " + pose_topic_name +"\33[0m")
        print("\33[94mpose_topic_type: " + pose_topic_type +"\33[0m")
        print("\33[94mcmd_vel_topic_name: " + cmd_vel_topic_name +"\33[0m")
        print("\33[94mpath_topic_name: " + path_topic_name +"\33[0m")
        print ("")
    except:
        print ("\33[41mProblem occurred when trying to read the parameters!: vec_field_alpha.py\33[0m")

    return

# ----------  ----------  ----------  ----------  ----------


# Main function
def vector_field():
    global freq
    global x_n, y_n, theta_n
    global pub_rviz_ref
    global has_path

    vel = Twist()

    i = 0

    # Init node
    rospy.init_node("vector_field")

    # Publisher for Twist message
    pub_cmd_vel = rospy.Publisher(cmd_vel_topic_name, Twist, queue_size=1)
    # Subscriber for pose
    if(pose_topic_type == "TFMessage"):
    	rospy.Subscriber(pose_topic_name, TFMessage, callback_tf)
    elif (pose_topic_type == "Pose"):
    	rospy.Subscriber(pose_topic_name, Pose, callback_pose)
    elif (pose_topic_type == "Odometry"):
    	rospy.Subscriber(pose_topic_name, Odometry, callback_odometry)
    elif (pose_topic_type == "PoseTurtle"):
        rospy.Subscriber(pose_topic_name, PoseTurtle, callback_turtlesim_pose)
    else :
    	print("\33[41mInvalid value for pose_topic_type!\33[0m")


    rospy.Subscriber(path_topic_name, Polygon, callback_path)


    # Publishers for rviz
    pub_rviz_ref = rospy.Publisher("/visualization_ref_vel", Marker, queue_size=1) #rviz marcador de velocidade de referencia

    rate = rospy.Rate(freq)

    print ("\33[92mVector field control started ...\33[0m\n")


    #Loop
    while not rospy.is_shutdown():

        # Count time
        i = i + 1
        time = i / float(freq)

        # Try to compute the vector field
        # try:
        if has_path:
            # Compute field
            [Vx_ref, Vy_ref] = vec_field_alpha(pos)

            # Compute a command of velocity
            [V_forward, w_z] = feedback_linearization(Vx_ref, Vy_ref)

            # Atribute values to the Twist message
            vel.linear.x = V_forward
            if(invert_motion_flag):
                vel.angular.z = -w_z
            else:
                vel.angular.z = w_z


            # Heuristic to detect that the end of the path was reached
            n_last = (len(path[0])+4)/3 - 1
            vec0 = [0.0, 0.0]
            vec0[0] = path[0][n_last]-pos[0]
            vec0[1] = path[1][n_last]-pos[1]
            vec = [0.0, 0.0]
            vec[0] = cos(rpy[2])*vec0[0] + sin(rpy[2])*vec0[1]
            vec[1] = -sin(rpy[2])*vec0[0] + cos(rpy[2])*vec0[1]
            number = vec[0]
            if(sqrt((pos[0]-path[0][n_last])**2+(pos[1]-path[1][n_last])**2) < 1.5):
                if(number < 0):

                    has_path = False
                    #Stop the robot
                    vel.linear.x = 0.0
                    vel.angular.z = 0.0

                    print ("\33[92mEnd of path reached\33[0m")

            # Publish velocity
            pub_cmd_vel.publish(vel)


            # Send markers to rviz
            # send_marker_to_rviz(pub_rviz_ref, Vx_ref, Vy_ref)

        # except:
        #     # This is due to the changes in the curve's change
        #     print "\33[93mProblem in the computation of the field !\33[0m"


        rate.sleep()


# ---------- !! ---------- !! ---------- !! ---------- !! ----------






# Funcao inicial
if __name__ == '__main__':

    # Frequency of field computation
    global freq
    freq = 20.0  # Hz

    # Robot position and orientation
    global pos, rpy
    pos = [0, 0, 0]
    rpy = [0, 0, 0]


    #Define parameters

    # Convergence intensity of the vector field
    global kf
    kf = 5.0

    # Constant relative to the feedback linearization controller
    global d
    d = 0.2

    # Parameter to invert the direction that the curve is being circulated
    global invert_direction
    invert_direction = False

    # Input parameters
    global vd

    # Flag to invert the sense of motion (move backwards)
    global invert_motion_flag
    invert_motion_flag = False

    # Names and type of topics
    global pose_topic_name, pose_topic_type, cmd_vel_topic_ame, path_topic_name
    pose_topic_name = "tf"
    pose_topic_type = "TFMessage"
    cmd_vel_topic_name = "cmd_vel"
    path_topic_name = "ref_path"


    # Sequence of points to be followed
    global path
    path = [[],[]]
    # Weights
    global W
    W = []
    # Start flag
    global has_path
    has_path = False

    #Load parameters
    load_parameters()

    try:
        vector_field()
    except rospy.ROSInterruptException:
        pass
