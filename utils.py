# Load standard modules
import numpy as np

# Load tudatpy modules
from tudatpy.interface import spice
from tudatpy import numerical_simulation
from tudatpy.numerical_simulation import environment_setup
from tudatpy.numerical_simulation import propagation_setup
from tudatpy.astro import element_conversion
from tudatpy.util import result2array

# Load pathlib
from pathlib import Path

####### FUNCTIONS #######

# -- retrieve_CR3BP_Transfer() --
# -> Retrieves the CR3BP transfer from the L2 Northern Halo to L1 Lyapunov orbit.
#    Computes respective patch points considering the number of total patch points selected.
def retrieve_CR3BP_Transfer(num_pp):
    # Basic checks for parameter compliance
    if not (num_pp>=40 and num_pp<=100): ValueError("Number of total patch points must be between 40 and 100")

    # Create file path to retrieve the CR3BP trajectory
    base = Path.cwd()
    file_path = base / "CR3BP_files" / "L2_L1_Transfer.txt"

    # Store full initial (CR3BP) solution
    q_nom = []
    t_nom = []
    with open(file_path, "r") as data:
        lines = data.readlines()
        for line in lines:
            x, y1, y2, y3, y4, y5, y6 = line.split()
            t_nom.append(float(x))
            q_nom.append(([float(y1), float(y2), float(y3), float(y4), float(y5), float(y6)]))

    # Convert to numpy arrays
    t_nom = np.array(t_nom)
    q_nom = np.array(q_nom)

    # Number of file lines
    N = len(t_nom)

    # Array of ids to define patch points
    idx = np.floor(np.linspace(0,N-1,num = num_pp+1)).astype(int)

    # Patch points and relative times of the first revolution
    q_nom_pp = q_nom[idx]
    t_nom_pp = t_nom[idx]

    # Return the results as numpy arrays
    return t_nom_pp, q_nom_pp, q_nom


# -- retrieve_CR3BP_PO() --
# -> Retrieves the CR3BP periodic orbit corresponding to the Lagrange point and orbit family selected.
#    Computes respective patch points considering the number of revolutions and patch points per
#    revolution selected.
def retrieve_CR3BP_PO(lagrange_point: str, orbit_family: str, num_rev: int, num_pp: int):
    # Basic checks for parameter compliance
    if not (lagrange_point == "L1" or lagrange_point == "L2"): ValueError("Lagrange Point must be either L1 or L2")
    if not (orbit_family == "NHalo" or orbit_family == "SHalo" or orbit_family == "Lyap" or orbit_family == "Vert"):
        ValueError("Orbit family must be either NHalo, SHalo, Lyap or Vert")
    if not num_rev>0: ValueError("Number of revolutions must be a positive integer")
    if not (num_pp>0 and num_pp<=10): ValueError("Number of patch points per revolution must be between 1 and 10")

    # Create file path to retrieve the CR3BP trajectory
    base = Path.cwd()
    file_name = lagrange_point + "_" + orbit_family + ".txt"
    file_path = base / "CR3BP_files" / file_name

    # Store full initial (CR3BP) solution
    q_nom = []
    t_nom = []
    with open(file_path, "r") as data:
        lines = data.readlines()
        for line in lines:
            x, y1, y2, y3, y4, y5, y6 = line.split()
            t_nom.append(float(x))
            q_nom.append(([float(y1), float(y2), float(y3), float(y4), float(y5), float(y6)]))

    # Convert to numpy arrays
    t_nom = np.array(t_nom)
    q_nom = np.array(q_nom)

    # Number of file lines
    N = len(t_nom)

    # Array of ids to define patch points
    idx = np.floor(np.linspace(0,N-1,num = num_pp+1)).astype(int)

    # Patch points and relative times of the first revolution
    q_pp_rev = q_nom[idx]
    t_pp_rev = t_nom[idx]

    # Create patch points and relative times for the number of revolutions selected
    t_nom_pp = []
    q_nom_pp = []
    t_offset = 0
    for i in range(num_rev):
        for q, t in zip(q_pp_rev[:-1], t_pp_rev[:-1]):
            q_nom_pp.append(q)
            t_nom_pp.append(t + t_offset)
        t_offset += t_pp_rev[-1]
    
    # Add final patch point
    q_nom_pp.append(q_pp_rev[-1])
    t_nom_pp.append(t_offset)

    q_nom_pp=np.array(q_nom_pp)
    t_nom_pp=np.array(t_nom_pp)

    # Return the results as numpy arrays
    return t_nom_pp, q_nom_pp, q_nom

# -- get_keplerian_elements_moon() --
# -> Retrieves the keplerian orbital elements of the Moon with respect to the Earth at the time instant provided
def get_keplerian_elements_moon(current_time, global_frame_orientation, body_settings):

    # Retrieve the cartesian state of the Moon at the time instant provided
    state_moon = spice.get_body_cartesian_state_at_epoch(
        target_body_name="Moon",
        observer_body_name="Earth",
        reference_frame_name=global_frame_orientation,
        aberration_corrections="NONE",
        ephemeris_time=current_time,
    )

    # Retrieve the gravitational parameter of the Earth
    mu_earth = body_settings.get("Earth").gravity_field_settings.gravitational_parameter

    # Convert the Moon's cartesian state into keplerian orbital elements with respect to the Earth
    moon_kepler = element_conversion.cartesian_to_keplerian(state_moon, mu_earth)

    # Return the Moon's keplerian orbital elements with respect to the Earth
    return moon_kepler

# -- get_Earth_Moon_distance() --
# -> Computes the magnitude of the distance between the Earth and Moon and the time instant provided
def get_Earth_Moon_distance(epochs, global_frame_origin, global_frame_orientation):

    l_EM = []

    if not isinstance(epochs, tuple):
        epochs = (epochs,)

    for epoch in epochs:
        # Retrieve the position of the Earth at the time instant provided
        pos_earth = spice.get_body_cartesian_position_at_epoch(
                target_body_name="Earth",
                observer_body_name=global_frame_origin,
                reference_frame_name=global_frame_orientation,
                aberration_corrections="NONE",
                ephemeris_time=epoch,
            )

        # Retrieve the position of the Moon at the time instant provided
        pos_moon= spice.get_body_cartesian_position_at_epoch(
                target_body_name="Moon",
                observer_body_name=global_frame_origin,
                reference_frame_name=global_frame_orientation,
                aberration_corrections="NONE",
                ephemeris_time=epoch,
            )
        
        # Compute the magnitude of the relative distance between the Earth and Moon
        l_EM.append(np.linalg.norm(pos_earth-pos_moon))
    
    l_EM = np.array(l_EM)
    return l_EM

# -- rotm_CR3BP_to_J2000() --
# -> Computes the rotation matrix from the CR3BP rotating frame to the J2000 inertial frame
def rotm_CR3BP_to_J2000(epochs, global_frame_origin, global_frame_orientation):

    R = []

    if not isinstance(epochs, tuple):
        epochs = (epochs,)

    for epoch in epochs:
        # Retrieve the initial Earth states
        state0_earth = spice.get_body_cartesian_state_at_epoch(
                target_body_name="Earth",
                observer_body_name=global_frame_origin,
                reference_frame_name=global_frame_orientation,
                aberration_corrections="NONE",
                ephemeris_time=epoch,
            )

        # Separate initial Earth states into position and velocity
        pos0_earth = state0_earth[0:3]
        vel0_earth = state0_earth[3:6]

        # Retrieve the initial Moon states
        state0_moon= spice.get_body_cartesian_state_at_epoch(
                target_body_name="Moon",
                observer_body_name=global_frame_origin,
                reference_frame_name=global_frame_orientation,
                aberration_corrections="NONE",
                ephemeris_time=epoch,
            )

        # Separate initial Moon states into position and velocity
        pos0_moon = state0_moon[0:3]
        vel0_moon = state0_moon[3:6]

        # Relative position and velocity vectors of the Moon with respect to the Earth are computed and normalized
        rel_pos0_moon = pos0_moon-pos0_earth
        r0 = rel_pos0_moon/np.linalg.norm(rel_pos0_moon)
        rel_vel0_moon = vel0_moon-vel0_earth
        v0 = rel_vel0_moon/np.linalg.norm(rel_vel0_moon)

        # The angular velocity of the Moon with respect to the Earth is computed and normalized
        rel_angvel0_moon = np.cross(r0,v0)
        w0 = rel_angvel0_moon/np.linalg.norm(rel_angvel0_moon)

        # Compute rotation matrix aligning X with r0, Z with w0, and Y completing the right-handed reference frame
        x = np.array([1.0, 0.0, 0.0])
        y = np.array([0.0, 1.0, 0.0])
        z = np.array([0.0, 0.0, 1.0])

        x_p = r0
        y_p = -np.cross(r0,w0)
        z_p = w0

        A = np.column_stack((x,y,z))
        B = np.column_stack((x_p,y_p,z_p))

        R.append(np.dot(B,A.T))

    R = np.array(R)

    # Return the rotation matrix array
    return R

# -- earth_moon_barycenter_state() --
# -> Computes the state of the Earth-Moon Barycenter (EMB) at the time instant provided
def earth_moon_barycenter_state(epoch, global_frame_origin, global_frame_orientation, body_settings):

    # Retrieve the Earth and Moon states at the time instant provided
    state_earth = spice.get_body_cartesian_state_at_epoch(
        target_body_name="Earth",
        observer_body_name=global_frame_origin,
        reference_frame_name=global_frame_orientation,
        aberration_corrections="NONE",
        ephemeris_time=epoch,
    )
    state_moon = spice.get_body_cartesian_state_at_epoch(
        target_body_name="Moon",
        observer_body_name=global_frame_origin,
        reference_frame_name=global_frame_orientation,
        aberration_corrections="NONE",
        ephemeris_time=epoch,
    )

    # Retrieve the gravitational parameters of the Earth and Moon
    mu_earth = body_settings.get("Earth").gravity_field_settings.gravitational_parameter
    mu_moon = body_settings.get("Moon").gravity_field_settings.gravitational_parameter

    # Compute barycenter state (position and velocity)
    barycenter_state = (mu_earth * state_earth + mu_moon * state_moon) / (mu_earth + mu_moon)

    # Return the state of the EMB
    return barycenter_state

def integrate_dynamics(q_i,tau_i,tau_f, central_bodies, acceleration_models, bodies_to_propagate, integrator_settings, bodies): # Takes dimensional quantities!

    # Create termination settings
    termination_condition = propagation_setup.propagator.time_termination(termination_time=tau_f,terminate_exactly_on_final_condition=True)

    # Create propagation settings
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        q_i,
        tau_i,
        integrator_settings,
        termination_condition,
    )

    # Create simulation object and propagate the dynamics
    dynamics_simulator = numerical_simulation.create_dynamics_simulator(
        bodies, propagator_settings
    )

    results = result2array(dynamics_simulator.propagation_results.state_history)
    t = results[:,0].copy()
    q = results[:,1:7].copy()
    
    return t,q # Outputs dimensional quantities!

def integrate_time(t,y,simulation_start_epoch,global_frame_origin, global_frame_orientation,mu_e,mu_m,T_SCALE,L_SCALE):
    y = y[0]
    L = get_Earth_Moon_distance(y*T_SCALE+simulation_start_epoch, global_frame_origin, global_frame_orientation)/L_SCALE
    dydt = np.sqrt((L**3)/(mu_e+mu_m))
    return dydt

def create_Q(pp_to_weight, weights, N_pp):
    Q = np.zeros((6*N_pp, 6*N_pp))

    for pp, weight in zip(pp_to_weight, weights):
        if pp>N_pp:
            ValueError("Patch-point to weight exceeds number of patch points")
        Q[6*pp:6*pp+3,6*pp:6*pp+3] = weight[0]*np.eye(3)
        Q[6*pp+3:6*pp+6, 6*pp+3:6*pp+6] = weight[1]*np.eye(3)

    return Q



######## FrameConverter class ########
# ~~ FrameConverter ~~~
# ~> Stores the necessary information to compute state conversions between the J2000 inertial frame 
# and the CR3BP rotating frame
class FrameConverter:

    # ~> Initialization requires the epoch for transformation and the non-dimensionalization scales
    # considered for the CR3BP (length, time, velocity, and acceleration)
    def __init__(self, epoch: float, mu_e: float, mu_m: float, global_frame_origin: str, 
                 global_frame_orientation: str, body_settings : environment_setup.BodySettings,
                 L_SCALE: float, T_SCALE: float, V_SCALE: float):
        
        # Step size used in the finite-differences
        step = 100

        # Epochs for finite-differences
        epochs = (epoch, epoch+step/2, epoch-step/2)

        d_array = get_Earth_Moon_distance(epochs, global_frame_origin, global_frame_orientation) 
        # Non-dimensional distance at the current time instant is retrieved.
        l = d_array[0]/ L_SCALE
        # Its (non-dimensional) derivatives are calculated via central finite-differences
        l_dot = (d_array[1]-d_array[2])/step / V_SCALE


        R_array = rotm_CR3BP_to_J2000(epochs, global_frame_origin, global_frame_orientation)
        # Rotation matrix between the two frames at the current time instant is retrieved.
        self.R = R_array[0]
        # Its (non-dimensional) derivatives are computed via central finite-differences 
        R_dot = (R_array[1]-R_array[2])/step * T_SCALE 

        # The state of the Earth-Moon Barycenter (EMB) at the current time instant is retrieved.
        q_emb = earth_moon_barycenter_state(epoch, global_frame_origin, global_frame_orientation, body_settings)
        # The position and velocity of the EMB are scaled accordingly
        self.r_emb = q_emb[0:3] / L_SCALE
        self.v_emb = q_emb[3:6] / V_SCALE

        # Auxiliary transformation coefficients
        self.m0 = l*self.R
        self.m00 = 1/l*self.R.T
        self.m1 = (l_dot*self.R + l*R_dot)
        self.tline = np.sqrt((mu_e+mu_m)/l**3)

    # -- CR3BP_to_J2000_inertial() --
    # -> Transforms the provided state from the CR3BP rotating frame into the J2000 inertial frame
    def CR3BP_to_J2000_inertial(self, q_E):
        
        # The J2000 state is initialized to zero
        q_J=np.zeros(len(q_E))

        # The position entries are transformed
        x_E = np.array(q_E[0:3])
        q_J[0:3] = self.r_emb + self.m0@x_E

        # If existent, the velocity entries are transformed
        if(len(q_E)>3):
            v_E = np.array(q_E[3:6])
            q_J[3:6] = self.v_emb + self.m1@x_E + self.m0*self.tline@v_E

        # The J2000 state is returned
        return q_J
    
    # -- J2000_inertial_to_CR3BP() --
    # -> Transforms the provided state from the J2000 inertial frame into the CR3BP rotating frame 
    def J2000_inertial_to_CR3BP(self, q_J):
        
        # The CR3BP state is initialized to zero
        q_E=np.zeros(len(q_J))

        # The position entries are transformed
        x_J = np.array(q_J[0:3])
        q_E[0:3] = self.m00 @(x_J - self.r_emb)

        # If existent, the velocity entries are transfromed
        if(len(q_E)>3):
            x_E = np.array(q_E[0:3])
            v_J = np.array(q_J[3:6])
            q_E[3:6] = 1/self.tline*self.m00 @(v_J - self.v_emb - self.m1@x_E)

        # The CR3BP state is returned
        return q_E   