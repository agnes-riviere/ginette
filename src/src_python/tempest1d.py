import numbers
import numpy as np
from scipy.linalg import expm, toeplitz
from filterpy.kalman import ExtendedKalmanFilter


class ModelProperties(object):
    """
    Properties of the model used to compute temperatures in 1D, including both
    mesh information and physical properties.
    
    Instantiated as ModelProperties(top_depth,
                                    bottom_depth,
                                    n_depths,
                                    Kt,
                                    Cw,
                                    Cs)
    Default values exist for every model property, so ModelProperties() returns
    an initialized object.

    Attributes
    ----------
    top_depth : float
        Depth of top boundary temperature, in meters
    bottom_depth : float
        Depth of bottom boundary temperature, in meters
    n_depths : int
        Number of depths at which temperatures are computed, including top and
        bottom boundaries
    Kt : float
        Thermal conductivity, in Watt/(meter * degree Celsius)
    Cw : float
        Heat capacity of water, in Joule/(meter^3 * degree Celsius)
    Cs : float
        Bulk heat capacity of water-sediment mixture, in Joule/(meter^3 *
        degree Celsius)
    nz : int
        Number of interior temperatures
    dz : float
        Mesh spacing
    zmin : float
        Depth of topmost interior temperature, in meters
    zmax : float
        Depth of bottommost interior temperature, in meters

    Notes
    -----
    The ModelProperties class holds all the properties needed to define the 1D
    problem. Many functions in tempest1d are written such that an instance of
    ModelProperties can be unpacked and used in the function arguments.
    
    Examples
    --------
    >>> import numpy as np
    >>> from tempest1d import ModelProperties
    >>> mp = ModelProperties(n_depths=101)
    >>> print(vars(mp))
    {'top_depth': 0, 'bottom_depth': 1, 'n_depths': 101, 'Kt': 2, 'Cw': 4182000.0, 'Cs': 2000000.0}
    """

    def __init__(self, top_depth=0, bottom_depth=1, n_depths=51, Kt=2,
                 Cw=4.182e6, Cs=2e6):
        self.top_depth = top_depth
        self.bottom_depth = bottom_depth
        self.n_depths = n_depths
        self.Kt = Kt
        self.Cw = Cw
        self.Cs = Cs
        # self.nz = n_depths - 2
        # assume uniformly spaced depths
        # self.dz = (bottom_depth - top_depth)/(n_depths - 1)
        # top internal temperature depth
        # self.zmin = top_depth + self.dz
        # bottom internal temperature depth
        # self.zmax = bottom_depth - self.dz

    @property
    def nz(self):
        """Number of interior temperatures"""
        return self.n_depths - 2

    @property
    def dz(self):
        """Mesh spacing, assuming uniformly spaced depths"""
        return (self.bottom_depth - self.top_depth)/(self.n_depths - 1)

    @property
    def zmin(self):
        """Depth of topmost interior temperature"""
        return self.top_depth + self.dz

    @property
    def zmax(self):
        """Depth of bottommost interior temperature"""
        return self.bottom_depth - self.dz

    @property
    def depths(self):
        """Array of all depth values, including boundaries"""
        return np.linspace(self.top_depth, self.bottom_depth, self.n_depths)


def _conduction_coefficient(mp):
    """
    Compute coefficient for conduction term, used in LHS matrices

    Parameters
    ----------
    mp : ModelProperties
        Model mesh and physical properties

    Returns
    -------
    float
        Conduction coefficient, in 1/second

    See Also
    --------
    _advection_coefficient
    """
    cc = mp.Kt/mp.Cs/(mp.dz**2)
    return cc


def _advection_coefficient(mp, q):
    """
    Compute coefficient for advection term, used in LHS matrices

    Parameters
    ----------
    mp : ModelProperties
        Model mesh and physical properties
    q : float
        Vertical specific discharge, in meters/second.
        Positive q is downward flux.

    Returns
    -------
    float
        Advection coefficient, in 1/second

    See Also
    --------
    _conduction_coefficient
    """
    ac = q*mp.Cw/mp.Cs/mp.dz
    return ac


def _A_expanded_old(mp, q):
    """
    Form matrix to compute dT/dt from T values including BCs
    dT/dt = _A_expanded * [T_top, T[:], T_bottom]

    Parameters
    ----------
    mp : ModelProperties
        Model mesh and physical properties
    q : float
        Vertical specific discharge, in meters/second.
        Positive q is downward flux.

    Returns
    -------
    A : nz by (nz + 2) sparse.dia.dia_matrix
        System dynamics matrix

    Notes
    -----
    This places boundary temperatures at -dz and zmax+dz
    """
    from scipy.sparse import diags
    diagonal = np.zeros(mp.nz)
    lower = np.zeros(mp.nz)
    upper = np.zeros(mp.nz)
    # conduction coefficient
    cc = _conduction_coefficient(mp)
    # advection coefficient
    ac = _advection_coefficient(mp, q)

    diagonal[:] = -2*cc-ac
    upper[:] = cc
    lower[:] = cc+ac
    A = diags(diagonals=[diagonal, upper, lower],
              offsets=[1, 2, 0], shape=(mp.nz, mp.nz + 2))
    return A


def _A_expanded(mp, q):
    """
    Form matrix to compute dT/dt from T values including BCs
    dT/dt = _A_expanded * [T_top, T[:], T_bottom]

    Parameters
    ----------
    mp : ModelProperties
        Model mesh and physical properties
    q : float
        Vertical specific discharge, in meters/second.
        Positive q is downward flux.

    Returns
    -------
    A : nz by (nz + 2) numpy ndarray
        System dynamics matrix

    Notes
    -----
    This places boundary temperatures at -dz and zmax+dz
    """
    # conduction coefficient
    cc = _conduction_coefficient(mp)
    # advection coefficient
    ac = _advection_coefficient(mp, q)

    # diagonal = np.zeros(mp.nz)
    # lower = np.zeros(mp.nz)
    # upper = np.zeros(mp.nz)
    # diagonal[:] = -2*cc-ac
    # upper[:] = cc
    # lower[:] = cc + ac
    # A = diags(diagonals=[diagonal, upper, lower],
    #           offsets=[1, 2, 0], shape=(mp.nz, mp.nz + 2))

    first_column = np.zeros(mp.nz)
    first_column[0] = cc + ac
    first_row = np.zeros(mp.nz + 2)
    first_row[0:3] = np.r_[cc + ac, -2*cc - ac, cc]
    A = toeplitz(first_column, first_row)
    return A


def _A(mp, q):
    """
    Form square matrix to compute dx/dt from x values
    _A @ x = dx/dt
    x = [T_top (BC), T1, ..., TN-1, T_bot (BC), q]

    Parameters
    ----------
    mp : ModelProperties
        Model mesh and physical properties
    q : float
        Vertical specific discharge, in meters/second.
        Positive q is downward flux.

    Returns
    -------
    A : (nz + 3) by (nz + 3) array
        Matrix that maps from x to dx/dt

    Examples
    --------
    Form the matrix:
    >>> mp = ModelProperties()
    >>> q = 1/24/60/60
    >>> A = _A(mp, q)
    >>> print(A.shape)
    (52, 52)
    """
    # A has all internal temps, 2 BCs, and q
    A = np.zeros((mp.nz+3, mp.nz+3))
    A[1:mp.nz+1, :mp.nz+2] = _A_expanded(mp, q)
    return A


def _expm(Adt):
    """
    Return exponential matrices F and B for _A, assuming constant boundary
    temperatures within a time step.
    
    This function serves as a helper, forming the exponential matrix and
    handling the indexing needed to return F and B.

    Parameters
    ----------
    Adt : ndarray
        A 2D square numpy array that computes dx from x values, as
        dx = Adt @ x. This matrix can be formed by multiplying the output of
        _A by a time step dt.

    Returns
    -------
    F : ndarray
        State transition matrix F
    B : ndarray
        Control function B

    See also
    --------
    _expm_interp
    """

    nz = Adt.shape[0] - 3
    # Most rows correspond to state variables
    # Two are for boundary temperatures, though
    x_rows = np.r_[np.arange(nz) + 1, -1]
    # same with columns
    x_cols = np.r_[np.arange(nz) + 1, -1]
    B_cols = np.r_[0, nz + 1]
    # Create 2D index tuples
    F_i = np.ix_(x_rows, x_cols)
    B_i = np.ix_(x_rows, B_cols)
    # Form exponential matrix
    expA = expm(Adt)
    # pick out non-boundary rows and columns
    F = expA[F_i]
    B = expA[B_i]
    return F, B


def _expm_interp(Adt):
    """
    Return exponential matrices F, B0, and B1 for _A,
    interpolating boundary temperatures across a time step.
    
    This function serves as a helper, forming the exponential matrix and
    handling the indexing needed to return F and B.

    Parameters
    ----------
    Adt : ndarray
        A 2D square numpy array that computes dx from x values, as
        dx = Adt @ x. This matrix can be formed by multiplying the output of
        _A by a time step dt.

    Returns
    -------
    A tuple containing
        F : ndarray
            State transition matrix F
        B0 : ndarray
             Control function corresponding to boundary temperatures at start
             of time step
        B1 : ndarray
             Control function corresponding to boundary temperatures at end of
             time step

    See also
    --------
    _expm
    """

    nz = Adt.shape[0] - 3
    nx = Adt.shape[0] - 2
    nu = 2
    I_u = np.zeros((nx+nu, nu))
    I_u[0, 0] = 1
    I_u[nz + 1, 1] = 1
    M = np.vstack([np.hstack([Adt, I_u]),
                   np.zeros((nu, nx + nu + nu))])
    # Form exponential matrix
    expM = expm(M)
    # Most rows correspond to state variables
    # Two are for boundary temperatures, though
    x_rows = np.r_[np.arange(nz) + 1, np.arange(nz + nu, nx + nu)]
    # same with columns
    x_cols = np.r_[np.arange(nz) + 1, np.arange(nz + nu, nx + nu)]
    B01_cols = np.r_[0, nz + 1]
    B1_cols = np.arange(nx + nu, nx + nu + nu)
    # Create 2D index tuples
    F_i = np.ix_(x_rows, x_cols)
    B01_i = np.ix_(x_rows, B01_cols)
    B1_i = np.ix_(x_rows, B1_cols)
    # pick out non-boundary rows and columns
    F = expM[F_i]
    B1 = expM[B1_i]
    B01 = expM[B01_i]
    B0 = B01 - B1
    return F, B0, B1


def _F_B(dt, mp, q, interp=True):
    """
    Form state transition matrices F and B

    Parameters
    ----------
    dt : float
        Size of the time step in seconds.
    mp : ModelProperties
        Model mesh and physical properties
    q : float
        Vertical specific discharge, in meters/second.
        Positive q is downward flux.
    interp : bool
        (Default value = True)
        If True, interpolate boundary temperatures linearly across time step.
        If False, use Tbc as constant values for boundary temperatures during
        time step

    Returns
    -------
    tuple
        If interp is True, returns (F, B0, B1) where F is the state transition
        matrix, B0 is the control function for boundary temperatures at the
        start of the time step, and B1 is the control function corresponding to
        boundary temperatures at the end of the time step. If interp is False,
        returns (F, B), where B is the control function for the constant
        boundary temperatures.
    """
    A = _A(mp, q)
    if interp:
        return _expm_interp(A*dt)
    else:
        return _expm(A*dt)


def simulate(x0, Tbc, dt, mp, q=None, interp=True, Tbc0=None):
    """
    Predict next state. State includes internal (non-boundary) temperatures and
    vertical discharge q. The simulation is done using a continuous time
    formulation.

    Parameters
    ----------
    x0 : length nz+1 ndarray
        Initial state. The array is comprised of internal temperatures and
        vertical discharge q. The temperatures are ordered from top to bottom,
        and the last element is q.
    Tbc : length 2 ndarray
        Boundary temperatures at the end of the time step.
    dt : float
        Size of the time step in seconds.
    mp : ModelProperties
        Model mesh and physical properties
    q : float
        Vertical specific discharge, in meters/second.
        Positive q is downward flux.
        If q is None, use the last element of x0 as next q. Usually q=None is a
        good choice. (Default value = None)
    interp : bool
        (Default value = True)
        If True, interpolate boundary temperatures linearly across time step.
        If False, use Tbc as constant values for boundary temperatures during
        time step
    Tbc0 : length 2 ndarray
        (Default value = None)
        Initial values for boundary temperatures. Only used if interp=True.

    Returns
    -------
    ndarray
        A length nz+1 array containing the computed state at the final time.
        The first nz elements are the internal temperatures, and the last
        element is the vertical discharge.

    Notes
    -----
    This places boundary temperatures at -dz and zmax+dz, where dz=zmax/nz.

    Examples
    --------
    Compute temperatures in the top meter of the subsurface after 10 minutes
    with a downward discharge velocity of 0.1 m/day, constant boundary
    temperatures of 20 deg C at the surface and 10 deg C at 1 m depth, and
    initial temperatures that vary linearly with depth.
    >>> import numpy as np
    >>> mp = ModelProperties()
    >>> q = 0.1/24/60/60
    >>> Tbc = np.array([20, 10])
    >>> dt = 10*60
    >>> initial_temperatures = np.linspace(Tbc[0], Tbc[1], mp.nz+2)
    >>> x0 = np.r_[initial_temperatures[1:-1], q]
    >>> xt = simulate(x0, Tbc, dt, mp, interp=False)
    >>> print(xt.shape)
    (50,)
    """
    if q is None:
        q = x0[-1]
    if interp:
        F, B0, B1 = _F_B(dt, mp, q, interp=interp)
        xbar = np.dot(F, x0) + np.dot(B0, Tbc0) + np.dot(B1, Tbc)
    else:
        F, B = _F_B(dt, mp, q, interp=interp)
        xbar = np.dot(F, x0) + np.dot(B, Tbc)
    return xbar


def _F_cd(F_linear, x0, Tbc, dt, mp, q=None, h=1e-8, interp=True, Tbc0=None):
    """
    Include dependency of temperature on discharge in Jacobian, dT_i/dq.
    Approximate derivative by central finite difference, as
    dT_i/dq ~= (T_i(q+h) - T_i(q-h))/2h

    Parameters
    ----------
    F_linear : 2D array
        The matrix exponential of A with dT/dq=0
            x_k = F_linear * x_k-1 + B_k * u_k
            A = [[0    , 0                  , 0    , 0],
                 [B_top, ...                , 0    , 0],
                 [0    , A_compact(**kwargs), 0    , 0],
                 [0    , ...                , B_bot, 0],
                 [0    , 0                  , 0    , 0],
                 [0    , 0                  , 0    , 0],
        F_linear = expm(A*dt)[1:nz+1 and nz+nu:, 1:nz+1 and nz+nu:]
    x0 : 1D array
        The current state vector
            x = [T_top, T_1, ..., T_nz, T_bot, q].T
    Tbc : length 2 ndarray
        Boundary temperatures at the end of the time step.
    dt : float
        Time step for Kalman filter
    mp : ModelProperties
        Model mesh and physical properties
    q : float
        (Default value = None)
        Vertical specific discharge, in meters/second.
        Positive q is downward flux.
        If q is None, use the last element of x0 as q
    h : float
        (Default value = 1e-8)
        The perturbation to apply to q for finite differencing
    interp : bool
        (Default value = True)
        If True, interpolate boundary temperatures linearly across time step.
        If False, use Tbc as constant values for boundary temperatures during
        time step
    Tbc0 : length 2 ndarray
        (Default value = None)
        Initial values for boundary temperatures. Only used if interp=True.

    Returns
    -------
    2D array
        The state transition matrix with dT/dq included.
    """
    if q is None:
        q = x0[-1]
    F = F_linear.copy()
    iq = mp.nz
    x0p = x0.copy()
    x0p[iq] += h
    xp = simulate(x0p, Tbc, dt, mp, q+h, interp=interp, Tbc0=Tbc0)
    x0m = x0.copy()
    x0m[iq] -= h
    xm = simulate(x0m, Tbc, dt, mp, q-h, interp=interp, Tbc0=Tbc0)
    F[:, iq] = (xp - xm)/(2*h)
    return F


class EKF(ExtendedKalmanFilter):
    """
    Extended Kalman Filter for estimating vertical discharge.
        State includes temperatures and discharge.
        Controls are top and bottom boundary temperatures.

    Parameters
    ----------
    dt : float
        time step for Kalman filter
    mp : ModelProperties
        Model mesh and physical properties
    interp : boolean
        (Default value = True)
        If True, interpolate boundary temperatures linearly across time step.
        If False, use Tbc as constant values for boundary temperatures during
        time step
    Tbc0 : length 2 ndarray
        (Default value = None)
        Initial values for boundary temperatures. Only used if interp=True.
    control_covariance : 2 by 2 array
        Covariance matrix for control inputs (temperatures)
    H : n_measurements by n_state array
        Matrix to interpolate state temperatures to measurement locations
    Hx : function
        Function to interpolate state temperatures to measurement locations
    HJacobian : n_measurements by n_state array
        Matrix to interpolate state temperatures to measurement locations
        Since the measurement function is linear, the Jacobian equals the matrix H.


    Examples
    --------
    >>> import numpy as np
    >>> from tempest1d import EKF, ModelProperties
    >>> mp = ModelProperties()
    >>> dt = 60*60
    >>> ekf = EKF(np.array([.2, .4, .7]), dt, mp, False)
    >>> Tbc = np.array([10, 20])
    >>> # initialize state
    >>> initial_temperatures = np.linspace(Tbc[0], Tbc[1], mp.nz+2)
    >>> q = 0.1/24/60/60
    >>> x0 = np.r_[initial_temperatures[1:-1], q]
    >>> ekf.x = x0
    >>> # initialize state covariance
    >>> ekf.P = np.eye(len(x0))
    >>> # process covariance
    >>> ekf.Q = np.eye(len(x0))*1e-4
    >>> ekf.Q[-1, -1] = 1e-14
    >>> # control covariance
    >>> measurement_std = 0.05
    >>> ekf.control_covariance = np.eye(2)*measurement_std**2
    >>> # measurement covariance
    >>> ekf.R = np.eye(3)*measurement_std**2
    >>> # predict next state
    >>> ekf.predict(Tbc)
    >>> # update state
    >>> ekf.update(np.r_[11.97, 13.95, 16.92])
    >>> # discharge estimate
    >>> print(ekf.x[-1])
    7.321506958157875e-07
    >>> # temperature estimates
    >>> print(ekf.Hx(ekf.x))
    [11.96908863 13.94999578 16.92090833]
    """
    def __init__(self, measure_depths, dt, mp, interp=True,
                 Tbc0=None, control_covariance=None):
        """
        Parameters
        ----------
        measure_depths : array
            Depths at which temperature measurements are taken
        dt : float
            time step for Kalman filter
        mp : ModelProperties
            Model mesh and physical properties
        interp : boolean
            (Default value = True)
            If True, interpolate boundary temperatures linearly across time step.
            If False, use Tbc as constant values for boundary temperatures during
            time step
        Tbc0 : length 2 ndarray
            (Default value = None)
            Initial values for boundary temperatures. Only used if interp=True.
        control_covariance : 2 by 2 array
            Covariance matrix for control inputs (temperatures)
        """
        n_state = mp.nz + 1
        self.dt = dt
        self.mp = mp
        n_measurements = len(measure_depths)
        super().__init__(dim_x=n_state, dim_z=n_measurements, dim_u=2)
        self.predict_x = self.predict_x
        self.interp = interp
        self.Tbc0 = Tbc0
        self.control_covariance = control_covariance
        z = mp.depths[1:-1]
        # Create interpolation matrix for observed temperatures
        H = np.zeros((n_measurements, n_state))
        for i_measurement, depth in enumerate(measure_depths):
            i_z = np.searchsorted(z, depth)
            if i_z == 0:
                H[i_measurement, i_z] = 1
            elif i_z >= n_state - 1:
                H[i_measurement, n_state - 1] = 1
            else:
                w = (depth - z[i_z-1])/(z[i_z] - z[i_z-1])
                H[i_measurement, i_z - 1] = 1 - w
                H[i_measurement, i_z] = w
        self.H = H
        self.HJacobian = lambda xx: self.H
        self.Hx = lambda xx: np.dot(self.H, xx)

    def predict_x(self, Tbc):
        """
        Predict next state with covariance using physics (F)
        Stores the next state prediction in member x

        Parameters
        ----------
        Tbc : ndarray
            a 2 element array of boundary temperatures
        """
        iq = self.mp.nz
        # get q from state
        q = self.x[iq]
        if self.interp:
            F_forward, Bk, Bk1 = _F_B(self.dt, self.mp, q, interp=self.interp)
            x_bar = (np.dot(F_forward, self.x) + np.dot(Bk, self.Tbc0) +
                     np.dot(Bk1, Tbc))
        else:
            F_forward, Bk = _F_B(self.dt, self.mp, q, interp=self.interp)
            x_bar = np.dot(F_forward, self.x) + np.dot(Bk, Tbc)
        # The F matrix for computing x_k+1 is different than the Jacobian
        # The difference is that the dependence of temperature on discharge is
        # only included in the Jacobian; it is zero in F_forward
        self.F = self._F_cd(F_forward, Tbc)
        self.x = x_bar
        self.Tbc0 = Tbc

    def predict(self, Tbc):
        """
        Prediction step
        Predict next state with covariance using physics (F)
        Stores the next state prediction in member x
        Updates state covariance in member P

        Parameters
        ----------
        Tbc : ndarray
            a 2 element array of boundary temperatures
        """
        iq = self.mp.nz
        # get q from state
        q = self.x[iq]
        if self.interp:
            F_forward, Bk, Bk1 = _F_B(self.dt, self.mp, q, interp=self.interp)
            x_bar = (np.dot(F_forward, self.x) + np.dot(Bk, self.Tbc0) +
                     np.dot(Bk1, Tbc))
        else:
            F_forward, Bk = _F_B(self.dt, self.mp, q, interp=self.interp)
            x_bar = np.dot(F_forward, self.x) + np.dot(Bk, Tbc)
        # The F matrix for computing x_k+1 is different than the Jacobian
        # The difference is that the dependence of temperature on discharge is
        # only included in the Jacobian; it is zero in F_forward
        self.F = self._F_cd(F_forward, Tbc)
        self.x = x_bar
        self.Tbc0 = Tbc

        self.P = np.dot(self.F, self.P).dot(self.F.T) + self.Q
        if self.control_covariance is not None:
            # Account for noise in the control input (Tbc)
            if self.interp:
                self.P += (np.dot(Bk, self.control_covariance).dot(Bk.T) +
                           np.dot(Bk1, self.control_covariance).dot(Bk1.T))
            else:
                self.P += np.dot(Bk, self.control_covariance).dot(Bk.T)

        # save prior
        self.x_prior = np.copy(self.x)
        self.P_prior = np.copy(self.P)

    def update(self, z, **kwargs):
        """
        Simplified update call. The measurement matrix H is linear, and it
        doesn't change between time steps. This function calls the internal
        update with the correct arguments for HJacobian and Hx.
        """
        super().update(z, self.HJacobian, self.Hx, **kwargs)

    def _F_cd(self, F_linear, Tbc, h=1e-8):
        """
        Include dependency of temperature on discharge in Jacobian, dT_i/dq.
        Approximate derivative by central finite difference, as
        dT_i/dq ~= (T_i(q+h) - T_i(q-h))/2h

        Parameters
        ----------
        F_linear : 2D array
            The matrix exponential of A with dT/dq=0
                x_k = F_linear * x_k-1 + B_k * u_k
                A = [[0    , 0                  , 0    , 0],
                     [B_top, ...                , 0    , 0],
                     [0    , A_compact(**kwargs), 0    , 0],
                     [0    , ...                , B_bot, 0],
                     [0    , 0                  , 0    , 0],
                     [0    , 0                  , 0    , 0],
            F_linear = expm(A*dt)[1:nz+1 and nz+nu:, 1:nz+1 and nz+nu:]
        Tbc : length 2 ndarray
            Boundary temperatures at the end of the time step.
        h : float
            (Default value = 1e-8)
            The perturbation to apply to q for finite differencing

        Returns
        -------
        2D array
            The state transition matrix with dT/dq included.
        """
        q = self.x[-1]
        F = F_linear.copy()
        iq = self.mp.nz
        x0p = self.x.copy()
        x0p[iq] += h
        xp = simulate(x0p, Tbc, self.dt, self.mp, q+h, interp=self.interp,
                      Tbc0=self.Tbc0)
        x0m = self.x.copy()
        x0m[iq] -= h
        xm = simulate(x0m, Tbc, self.dt, self.mp, q-h, interp=self.interp,
                      Tbc0=self.Tbc0)
        F[:, iq] = (xp - xm)/(2*h)
        return F


def run_EKF(ekf, measurements, T_top, T_bottom, 
            dt=None,
            Qc=None,
            return_full_P=False):
    """
    Perform Kalman filter estimates for measurements at multiple times.

    Parameters
    ----------
    ekf : EKF object
        Initialized EKF object
    measurements : nz by nt array
        Temperature measurements at all times and depths.
    T_top : 1D array 
        Boundary temperatures at top of model.
    T_bottom : 1D array 
        Boundary temperatures at bottom of model.
    dt : float or 1D array
        (Default value = None)
        Lengths of time intervals between measurements for Kalman filter.
        If a float, every time interval is assumed to be the same.
        If None, ekf.dt is used.
    Qc : nx by nx array
        (Default value = None)
        Process covariance before integrating over time.
        If not None, Q = Qc*dt
        If None, Q = ekf.Q
    return_full_P : boolean
        (Default value = False)
        Whether to return the full state covariance matrix for all time steps,
        or only the discharge variances. The full state covariance matrix is
        needed for RTS smoothing.

    Returns
    -------
    Tuple of (x_ekf, y_ekf, P_ekf)
        x_ekf : nt by nx array
            Posterior state estimates (after update) for all time steps.
        y_ekf : nt by nz array
            Measurement residuals for all time steps.
        P_ekf : nt [by nx by nx] array
            If return_full_P, posterior covariance for all time steps.
            If not, posterior discharge variance for all time steps.
    
    Examples
    --------
    >>> import numpy as np
    >>> from tempest1d import EKF, ModelProperties, run_EKF
    >>> mp = ModelProperties()
    >>> # One temperature profile per hour
    >>> dt = 60*60
    >>> # Three depths: 0.2 m, 0.4 m, and 0.7 m
    >>> ekf = EKF(np.array([.2, .4, .7]), dt, mp, False)
    >>> # boundary temperatures through four time steps
    >>> T_top = np.array([17.9, 18.1, 18.7, 19.3])
    >>> T_bottom = np.array([15.6, 15.6, 15.6, 15.6])
    >>> # estimate initial state by interpolating boundary temperatures
    >>> initial_temperatures = np.linspace(T_top[0], T_bottom[0], mp.nz+2)
    >>> q = 0.1/24/60/60
    >>> x0 = np.r_[initial_temperatures[1:-1], q]
    >>> ekf.x = x0
    >>> # initialize state covariance
    >>> ekf.P = np.eye(len(x0))
    >>> # process covariance
    >>> ekf.Q = np.eye(len(x0))*1e-4
    >>> ekf.Q[-1, -1] = 1e-14
    >>> # control covariance
    >>> measurement_std = 0.05
    >>> ekf.control_covariance = np.eye(2)*measurement_std**2
    >>> # measurement covariance
    >>> ekf.R = np.eye(3)*measurement_std**2
    >>> # measurements at each depth
    >>> measurements = np.array([[17.7, 17.7, 17.8, 17.9],
    >>>                          [17.4, 17.4, 17.4, 17.5],
    >>>                          [16.8, 16.8, 16.8, 16.8]])
    >>> x_ekf, y_ekf, P_ekf = run_EKF(ekf, measurements, T_top, T_bottom, 
    >>>     dt=dt,
    >>>     Qc=Q/dt,
    >>>     return_full_P=False):
    >>> # discharge estimate
    >>> print(ekf.x[-1])
    7.321506958157875e-07
    >>> # temperature estimates
    >>> print(ekf.Hx(ekf.x))
    [11.96908863 13.94999578 16.92090833]
    """
    nx = ekf.dim_x
    nz = ekf.dim_z
    z = measurements.reshape((nz, -1))
    nt = z.shape[1]

    # parse dt
    if isinstance(dt, numbers.Number):
        dt = [dt]*nt
        ekf.dt = dt[0]
    elif dt is None:
        dt = [ekf.dt]*nt
    elif len(dt) != nt-1:
        raise ValueError('Number of measurements must be one more that len(dt)')
    else:
        ekf.dt = dt[0]

    # TODO: parse Q

    # check for nans among inputs
    # control_nans = np.any(np.isnan([T_top, T_bottom]), axis=0)
    # measurement_nans = np.any(np.isnan(z), axis=0)

    # initialize results arrays
    # TODO: base output size on control_nans
    x_ekf = np.zeros((nt, nx))
    x_ekf[0, :] = ekf.x
    y_ekf = np.zeros((nt, nz))
    y_ekf[0, :] = z[:, 0] - ekf.Hx(ekf.x)
    if return_full_P:
        P_ekf = np.zeros((nt, nx, nx))
        P_ekf[0] = ekf.P
    else:
        # only return discharge variance
        P_ekf = np.zeros(nt)
        P_ekf[0] = ekf.P[-1, -1]

    # run Kalman filter at each time step
    for i_t in range(1, nt):
        ekf.dt = dt[i_t - 1]
        if Qc is not None:
            ekf.Q = Qc*ekf.dt
        # if control_nans[i_t]:
            # ekf.dt = ekf.dt + dt[i_t - 1]
        # else:
        # predict
        ekf.predict(np.r_[T_top[i_t], T_bottom[i_t]])
        # if not measurement_nans[i_t]:
        # update
        z_t = z[:, i_t]
        ekf.update(z_t)
        # if i_t < nt - 1:
            # ekf.dt = dt[i_t]
        # save to results arrays
        x_ekf[i_t, :] = ekf.x
        y_ekf[i_t] = ekf.y
        if return_full_P:
            P_ekf[i_t, :, :] = ekf.P
        else:
            P_ekf[i_t] = ekf.P[-1,-1]
    return (x_ekf, y_ekf, P_ekf)


def run_RTS(ekf, x_ekf, P_ekf, measurements, T_top, T_bottom, 
            dt=None,
            Qc=None,
            return_full_P=False):
    """
    Perform Rauch-Tung-Striebel smoothing for measurements at multiple times.

    Parameters
    ----------
    ekf : EKF object
        Initialized EKF object
    x_ekf : 1D array of nx elements
        Posterior state estimates obtained from, e.g., run_EKF()
    P_ekf : nt by nx by nx array
        Posterior state covariance for all time steps obtained from, e.g.,
        run_EKF()
    measurements : nz by nt array
        Temperature measurements at all times and depths.
    T_top : 1D array 
        Boundary temperatures at top of model.
    T_bottom : 1D array 
        Boundary temperatures at bottom of model.
    dt : float or 1D array
        (Default value = None)
        Lengths of time intervals between measurements for Kalman filter.
        If a float, every time interval is assumed to be the same.
        If None, ekf.dt is used.
    Qc : nx by nx array
        (Default value = None)
        Process covariance before integrating over time.
        If not None, Q = Qc*dt
        If None, Q = ekf.Q
    return_full_P : boolean
        (Default value = False)
        Whether to return the full state covariance matrix for all time steps,
        or only the discharge variances. The full state covariance matrix is
        needed for RTS smoothing.

    Returns
    -------
    Tuple of (x_rts, y_rts, P_rts)
        x_rts : nt by nx array
            Smoothed posterior state estimates (after update) for all time steps.
        y_rts : nt by nz array
            Measurement residuals for all time steps after smoothing.
        P_rts : nt [by nx by nx] array
            If return_full_P, posterior covariance for all time steps after
            smoothing.  If not, posterior discharge variance for all time steps
            after smoothing.
    """
    nx = ekf.dim_x
    nz = ekf.dim_z
    z = measurements.reshape((nz, -1))
    nt = z.shape[1]
    # initialize results arrays
    x_rts = np.zeros((nt, nx))
    y_rts = np.zeros((nt, nz))
    P_rts = np.zeros((nt, nx, nx))
    x_rts[-1, :] = x_ekf[-1, :]
    P_rts[-1, :, :] = P_ekf[-1, :, :]

    # parse dt
    if isinstance(dt, numbers.Number):
        dt = [dt]*nt
        ekf.dt = dt[0]
    elif dt is None:
        dt = [ekf.dt]*nt
    elif len(dt) != nt-1:
        raise ValueError('Number of measurements must be one more that len(dt)')
    else:
        ekf.dt = dt[0]

    # RTS smoothing, starting from most recent time step
    for i_t in range(nt-2, -1, -1):
        q_i = x_ekf[i_t, -1]
        dt_i = dt[i_t]
        if Qc is not None:
            ekf.Q = Qc*dt_i
        Tbc = np.r_[T_top[i_t + 1], T_bottom[i_t + 1]]
        Tbc0 = np.r_[T_top[i_t], T_bottom[i_t]]
        # Form Jacobian and other transition matrices
        F_linear, Bk, Bk1 = _F_B(dt_i, ekf.mp, q_i, interp=True)
        x_bar = (np.dot(F_linear, x_ekf[i_t]) + np.dot(Bk, Tbc0) +
                 np.dot(Bk1, Tbc))
        F = _F_cd(F_linear, x_ekf[i_t], 
                  Tbc,
                  dt_i,
                  ekf.mp,
                  q=q_i,
                  interp=ekf.interp,
                  Tbc0=Tbc0
                 )

        P_bar = np.dot(F, P_ekf[i_t]).dot(F.T) + ekf.Q
        P_bar += (np.dot(Bk, ekf.control_covariance).dot(Bk.T) +
                   np.dot(Bk1, ekf.control_covariance).dot(Bk1.T))
        G = np.dot(np.dot(P_ekf[i_t], F.T), 
                   np.linalg.inv(P_bar)
                  )
        x_rts[i_t, :] = x_ekf[i_t] + np.dot(G, x_rts[i_t+1] - x_bar)
        P_rts[i_t, :] = P_ekf[i_t] + np.dot(np.dot(G, P_rts[i_t+1] - P_bar), G.T)

    # Compute residuals
    z_rts = ekf.Hx(x_rts.T)
    y_rts = (measurements - z_rts).T
    if return_full_P:
        P_rts = P_rts
    else:
        P_rts = P_rts[:, -1, -1]
    return (x_rts, y_rts, P_rts)

