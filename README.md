# Vector Field Control
Control codes based on artificial vector fields

---------------------
This repository performs the control of a differential robot by using artificial vector fields. To use it, just clone it to your catkin_ws/src folder:

```bash
$ git clone https://github.com/adrianomcr/vector_field_control.git
```



## Vector field control (Alphas):

The node `vec_field_alpha.py` is used to control the robot to follow a path represented as the zero level set of a function alpha. It requires the global pose of the robot. The current alpha function implemented is:

![formula](https://render.githubusercontent.com/render/math?math=\alpha(x,y)=\left(\left(\frac{x-c_x}{a}\right)^\gamma%2B\left(\frac{y-c_y}{b}\right)^\gamma\right)^{\frac{1}{\gamma}}-1)


### Topics

- `pose_topic_name`  (message type: `tf2_msgs/TFMessage`): Subscribe to this topic to get the pose of the robot. To be selected with parameter `pose_topic_type`;
- `pose_topic_name`  (message type: `geometry_msgs/Pose`): Subscribe to this topic to get the pose of the robot. To be selected with parameter `pose_topic_type`;
- `pose_topic_name`  (message type: `nav_msgs/Odometry`): Subscribe to this topic to get the pose of the robot. To be selected with parameter `pose_topic_type`;
- `pose_topic_name`  (message type: `turtlesim_msgs/Pose`): Subscribe to this topic to get the pose of the robot. To be selected with parameter `pose_topic_type`;
- `cmd_vel_topic_name`  (message type: `geometry_msgs/Twist`): Publish at this topic a velocity command (forward velocity and an angular velocity)

Note: The type of topic whose name  is in the parameter `pose_topic_name` should be selected with the parameter `pose_topic_type`. This enables the code to get the pose from different message types.

### Parameters

The list of parameters is on the file inside the `config` folder in the file `control_params.yaml`. They are:

- `d_feedback` (`float`): Distance the control point is moved forward from the robot's center;
- `K_F` (`float`): Convergence gain of the vector field;
- `vd` (`float`): Reference forward speed for the espeleorobo;
- `a` (`float`): Stretching of the curve on x axis
- `b` (`float`): Stretching of the curve on y axis
- `cx` (`float`): Center of the curve on x axis
- `cy` (`float`): Center of the curve on y axis
- `gamma` (`int`): Parameter of the curve's shape. Should be an even integer
- `invert_direction` (`bool`): Flag to invert the sense of the curve's circulation
- `invert_motion_flag` (`bool`): Flag to invert the motion of the espeleorobo (move backwards);
- `pose_topic_name` (`string`): Name of the topic in which the pose will be obtained;
- `pose_topic_type` (`string`): Type of the topic in which the pose will be obtained (options: TFMessage, Pose, Odometry);
- `cmd_vel_topic_name` (`string`): Name of the topic in which the forward and angular velocities will be published.



### Usage

To use this node just add the following lines in your launch file:

```xml
<!-- Run the node that controls the robot with vector fields -->
<node pkg = "vector_field_control" name = "vector_field" type = "vec_field_alpha.py" args="" output="screen">
  <rosparam command="load" file="$(find vector_field_control)/config/control_params.yaml" />
</node>
```

You can simply use the file `vector_field.launch`:

```bash
$ roslaunch vector_field_control vector_field.launch
```

#### Note
 This node was developed to control the EspeleoRobô, a platform developed by Vale S.A., the biggest mining company in the world. The models necessary to simulate the EspeleoRobô are not yet available for the general public. Thus, the config parameters are set to run an example with the turtlesim simulator.






## Vector field control (Distance field):

The node `vec_field_dist.py` is used to control the robot to follow a path represented by a parametric equation. The field uses the Euclidean distance. It requires the global pose of the robot. The parametric rquation is provides with two strings.


### Topics

- `pose_topic_name`  (message type: `tf2_msgs/TFMessage`): Subscribe to this topic to get the pose of the robot. To be selected with parameter `pose_topic_type`;
- `pose_topic_name`  (message type: `geometry_msgs/Pose`): Subscribe to this topic to get the pose of the robot. To be selected with parameter `pose_topic_type`;
- `pose_topic_name`  (message type: `nav_msgs/Odometry`): Subscribe to this topic to get the pose of the robot. To be selected with parameter `pose_topic_type`;
- `pose_topic_name`  (message type: `turtlesim_msgs/Pose`): Subscribe to this topic to get the pose of the robot. To be selected with parameter `pose_topic_type`;
- `cmd_vel_topic_name`  (message type: `geometry_msgs/Twist`): Publish at this topic a velocity command (forward velocity and an angular velocity)

Note: The type of topic whose name  is in the parameter `pose_topic_name` should be selected with the parameter `pose_topic_type`. This enables the code to get the pose from different message types.

### Parameters

The list of parameters is on the file inside the `config` folder in the file `control_params_dist.yaml`. They are:

- `d_feedback` (`float`): Distance the control point is moved forward from the robot's center;
- `K_F` (`float`): Convergence gain of the vector field;
- `vd` (`float`): Reference forward speed for the espeleorobo;
- `N` (`int`): Number of samples of he curve
- `s1` (`float`): Minimum value of the paremeter s
- `s2` (`float`): Maximum value of the paremeter s
- `rx_str` (`string`): Parametric representation of the curve on x axis
- `ry_str` (`string`): Parametric representation of the curve on y axis
- `gamma` (`int`): Parameter of the curve's shape. Should be an even integer
- `invert_direction` (`bool`): Flag to invert the sense of the curve's circulation
- `invert_motion_flag` (`bool`): Flag to invert the motion of the espeleorobo (move backwards);
- `pose_topic_name` (`string`): Name of the topic in which the pose will be obtained;
- `pose_topic_type` (`string`): Type of the topic in which the pose will be obtained (options: TFMessage, Pose, Odometry);
- `cmd_vel_topic_name` (`string`): Name of the topic in which the forward and angular velocities will be published.



### Usage

To use the distance vector field just add the following lines in your launch file:

```xml
<!-- Run the node that controls the robot with vector fields (Distance field) -->
<node pkg = "vector_field_control" name = "vector_field" type = "vec_field_dist.py" args="" output="screen">
  <rosparam command="load" file="$(find vector_field_control)/config/control_params_dist.yaml" />
</node>
```

You can simply use the file `vector_field_dist.launch`:

```bash
$ roslaunch vector_field_control vector_field_dist.launch
```




## Exmples

Example of the vector field based on alphas:

```bash
$ roscore
$ rosrun turtlesim turtlesim_node
$ roslaunch vector_field_control vector_field.launch
```

If everything is ok you should see the following result:

![image](https://github.com/adrianomcr/vector_field_control/blob/master/images/turtle.png)



Example of the vector field based on the Euclidean distance:

```bash
$ roscore
$ rosrun turtlesim turtlesim_node
$ roslaunch vector_field_control vector_field_dist.launch
```

If everything is ok you should see the following result:

![image](https://github.com/adrianomcr/vector_field_control/blob/master/images/turtle_dist.png)

















## Contact:

Adriano Rezende

adrianomcr18@gmail.com
