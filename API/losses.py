import numpy as np
from scipy.interpolate import interp1d
import math


def Npr(Re, ks, width):
    coeff = np.array([[3.1, 3.1, 3.1, 3.1, 3.1, 3.1],
                      [2.4, 2.4, 2.4, 2.4, 2.4, 2.4],
                      [2.1, 1.8, 1.8, 1.8, 1.8, 1.8],
                      [2.1, 1.6, 1.6, 1.6, 1.6, 1.6],
                      [2.1, 1.35, 1.35, 1.35, 1.35, 1.35],
                      [2.1, 1.25, 1.05, 1.05, 1.05, 1.05],
                      [2.1, 1.25, 1, 0.9, 0.9, 0.9],
                      [2.1, 1.25, 1, 0.75, 0.75, 0.75],
                      [2.1, 1.25, 1, 0.7, 0.65, 0.65]])
    Re_coeff = np.array([1.0e4, 5.0e4, 1.0e5, 5.0e5, 1.0e6, 5.0e6])
    ks_coeff = np.array([10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02])

    f1 = [interp1d(Re_coeff, i, fill_value='extrapolate') for i in coeff]
    y1 = [func(Re) for func in f1]
    return interp1d(ks_coeff, y1, fill_value='extrapolate')(1000 * ks / width).tolist()


def Fl(inlet_angle, outlet_angle):
    outlet_range = np.arange(10, 100, 10)
    inlet_range = np.arange(40, 150, 10)
    points_num = np.array([4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4])
    Fl_coeff = np.array([[12.45, 13.7, 14.6, 15.1, -1, -1, -1, -1, -1],
                         [12.5, 13.3, 13.9, 13.95, 13.6, -1, -1, -1, -1],
                         [12.25, 13, 13.25, 12.9, 12.2, 11.15, -1, -1, -1],
                         [12.0, 13.0, 12.9, 12.15, 10.9, 9.4, 8.0, -1, -1],
                         [12, 12.75, 12.25, 11.2, 9.6, 7.7, 5.75, 3.75, -1],
                         [12, 12.5, 11.83, 10.275, 8.2, 6.05, 3.85, 1.8, 0],
                         [11.5, 11.15, 10.15, 8.7, 6.5, 4, 1.9, 0, -1],
                         [10.75, 10, 8.55, 6.6, 4.35, 2.15, 0, -1, -1],
                         [10, 8.55, 6.66, 4.6, 2.4, 0, -1, -1, -1],
                         [9, 7.2, 4.9, 2.35, 0, -1, -1, -1, -1],
                         [8, 5.5, 2.8, 0, -1, -1, -1, -1, -1]])
    f1 = [interp1d(outlet_range[:points_num[j]], i[:points_num[j]], fill_value='extrapolate') for j, i in
          enumerate(Fl_coeff)]
    y1 = [func(outlet_angle) for func in f1]
    return interp1d(inlet_range, y1, fill_value='extrapolate')(inlet_angle).tolist()


def Cr(inlet_angle, outlet_angle, pitch, width):
    Cr_coeff = np.array([[1.03, 1.02, 1.035, 1.05, 1.12, 1.21, 1.32, 1.53, 1.65, 1.95],
                         [1.05, 1.05, 1.07, 1.1, 1.15, 1.25, 1.5, 1.68, 2, 2.4],
                         [1.1, 1.105, 1.12, 1.14, 1.2, 1.34, 1.6, 1.9, 2.3, 2.8]])
    s_b_ratio_range = np.array([1.25, 0.8, 0.4])
    angles_ratio_range = np.arange(-0.2, 0.8, 0.1)
    s_b = pitch / width
    angle_ratio = math.sin(math.radians(outlet_angle)) / math.sin(math.radians(inlet_angle))
    f1 = [interp1d(angles_ratio_range, i, fill_value='extrapolate') for i in Cr_coeff]
    y1 = [func(angle_ratio) for func in f1]
    return interp1d(s_b_ratio_range, y1, fill_value='extrapolate')(s_b).tolist()


def Xpb(inlet_angle, outlet_angle, pitch, width):
    Fl_mod_range = np.arange(2, 12, 1)
    Cr_range = np.array([1.0, 1.1, 1.3, 1.5, 2.0, 5.0])
    points_num = np.array([8, 9, 9, 10, 10, 10])
    Xpb_coeff = np.array([[1.5, 1.5, 1.5, 1.55, 1.59, 1.78, 2.39, 3.8, -1, -1],
                          [1.2, 1.2, 1.2, 1.2, 1.2, 1.27, 1.48, 2.05, 2.9, -1],
                          [1, 1, 1, 1, 1, 1.06, 1.14, 1.39, 1.89, -1],
                          [0.88, 0.88, 0.88, 0.88, 0.88, 0.89, 0.95, 1.14, 1.53, 2.12],
                          [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.75, 0.86, 1.12, 1.5],
                          [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.45, 0.507, 0.725]])
    Fl_mod = Fl(inlet_angle, outlet_angle) * pitch / width
    Cr_val = Cr(inlet_angle, outlet_angle, pitch, width)
    f1 = [interp1d(Fl_mod_range[:points_num[j]], i[:points_num[j]], fill_value='extrapolate') for j, i in
          enumerate(Xpb_coeff)]
    y1 = [func(Fl_mod) for func in f1]
    return interp1d(Cr_range, y1, fill_value='extrapolate')(Cr_val).tolist()


def Npt(te_radius, pitch, outlet_angle):
    te_pitch_range = np.arange(0, 0.13, 0.02)
    outlet_angle_range = np.array([10, 15, 20, 30, 50])
    Npt_coeff = np.array([[1, 1.15, 1.34, 1.61, 2.02, 2.52, -1],
                          [1, 1.11, 1.225, 1.36, 1.53, 1.72, 1.96],
                          [1, 1.09, 1.16, 1.26, 1.36, 1.5, 1.6],
                          [1, 1.08, 1.15, 1.225, 1.3, 1.36, 1.46],
                          [1, 1.07, 1.135, 1.19, 1.25, 1.31, 1.375]])
    te_pitch = te_radius / pitch
    f1 = [interp1d(te_pitch_range[:-1], i[:-1], fill_value='extrapolate') if j == 0 else
          interp1d(te_pitch_range, i, fill_value='extrapolate') for j, i in enumerate(Npt_coeff)]
    y1 = [func(te_pitch) for func in f1]
    return interp1d(outlet_angle_range, y1, fill_value='extrapolate')(outlet_angle).tolist()


def dXpte(te_radius, pitch):
    te_pitch_range = np.arange(0, 0.13, 0.02)
    dXpte_coeff = np.array([0, 0.06, 0.17, 0.385, 0.6, 1.1, 1.85])
    return interp1d(te_pitch_range, dXpte_coeff, fill_value='extrapolate', kind='cubic')(te_radius / pitch).tolist()


def dMach():
    return 0.0


def dBBr():
    return 0.0


# figure 2.23
def CR_inc(inlet_angle, outlet_angle, pitch, width, throat):
    CR_range = np.arange(1, 3.1, 0.5)
    exit_angle_range = np.arange(60, 9, -10)
    CR_inc_coeff = np.array([[0, 2, 5, 13, 18.5],
                             [0, 2, 5, 13, 18.5],
                             [0, 2, 5, 13, 18.5],
                             [-1, 1, 5, 13, 18.5],
                             [-3.5, 0.5, 5, 13, 18.5],
                             [-3.5, 0.5, 5, 13, 18.5]])
    f1 = [interp1d(CR_range, i, fill_value='extrapolate') for i in CR_inc_coeff]
    CR_ratio = Cr(inlet_angle, outlet_angle, pitch, width)
    y1 = [func(CR_ratio) for func in f1]
    exit_angle_val = math.degrees(math.asin(throat / pitch))
    return interp1d(exit_angle_range, y1, fill_value='extrapolate')(exit_angle_val).tolist()


# figure 2.22
def PW_inc_positive(pitch, width, throat):
    PW_range = np.arange(0.3, 1.8, 0.2)
    exit_angle_range = np.arange(30, 75, 10)
    PW_inc_coeff = np.array([[20.5, 18.5, 7.5, -8, -18, -23, -25.5, -24],
                             [15.5, 14, 5, -5.5, -14.5, -20.5, -22.5, -20],
                             [10, 9.5, 4, -5, -11.5, -17.5, -19.5, -14.5],
                             [7, 6.5, 3, -4, -9, -14, -15, -10],
                             [4, 4, 2, -2.5, -6.5, -10.5, -11, -5]])
    f1 = [interp1d(PW_range, i, fill_value='extrapolate') for i in PW_inc_coeff]
    y1 = [func(pitch / width) for func in f1]
    exit_angle_val = math.degrees(math.asin(throat / pitch))
    return interp1d(exit_angle_range, y1, fill_value='extrapolate')(exit_angle_val).tolist()


# figure 2.25
def PW_inc_negative(pitch, width, throat):
    PW_range = np.arange(0.3, 1.6, 0.2)
    exit_angle_range = np.arange(10, 55, 10)
    PW_inc_coeff = np.array([[-17.5, -10, -3, 2.5, 5.5, 5.5, 5.5],
                             [-21, -11.5, -4, 2.5, 5.5, 5.5, 5.5],
                             [-25, -13, -4.5, 2.5, 5.5, 5.5, 5.5],
                             [-30, -16, -5, 2.5, 5.5, 5.5, 5.5],
                             [-35, -20, -6, 2.5, 5.5, 5.5, 5.5]])
    f1 = [interp1d(PW_range, i, fill_value='extrapolate') for i in PW_inc_coeff]
    y1 = [func(pitch / width) for func in f1]
    exit_angle_val = math.degrees(math.asin(throat / pitch))
    return interp1d(exit_angle_range, y1, fill_value='extrapolate')(exit_angle_val).tolist()


# figure 2.27
def FI_min(blade_inlet_angle, pitch, width):
    inlet_angle_range = np.arange(30, 111, 10)
    PW_range = np.array([1.1, 0.9, 0.7, 0.5, 0.4])
    FI_coeff = np.array([[1, 1, 1.04, 1.1, 1.2, 1.35, 1.45, 1.6, 1.77],
                         [0.8, 0.82, 0.88, 0.92, 1, 1.1, 1.25, 1.4, 1.55],
                         [0.59, 0.61, 0.65, 0.7, 0.8, 0.91, 1.05, 1.2, 1.32],
                         [0.35, 0.4, 0.44, 0.48, 0.56, 0.67, 0.8, 0.95, 1.1],
                         [0.21, 0.25, 0.3, 0.35, 0.43, 0.51, 0.61, 0.8, 0.95]])
    f1 = [interp1d(inlet_angle_range, i, fill_value='extrapolate') for i in FI_coeff]
    y1 = [func(blade_inlet_angle) for func in f1]
    return interp1d(PW_range, y1, fill_value='extrapolate')(pitch / width).tolist()


# figure 2.21
def i_val_inlet_lower_90(blade_inlet_angle, pitch, throat):
    blade_inlet_angle_range = np.arange(50, 125, 10)
    exit_angle_range = np.arange(60, 5, -10)
    i_angle_coeff = np.array([[5, 11.5, 12.5, 12.4, 12.4, 12.3, 12.1, 11.5],
                              [3, 9, 15, 16.5, 17, 16.5, 15, 12.5],
                              [1, 7.5, 14.5, 21, 25.5, 26, 24, 20],
                              [-1, 5, 11.5, 18.5, 25.5, 30, 30.5, 27.5],
                              [-5, 2, 8, 15, 20.5, 24.5, 27, 27.5],
                              [-8, -2, 2.5, 7, 11.5, 16, 19.5, 22]])
    f1 = [interp1d(blade_inlet_angle_range, i, fill_value='extrapolate') for i in i_angle_coeff]
    y1 = [func(blade_inlet_angle) for func in f1]
    return interp1d(exit_angle_range, y1, fill_value='extrapolate')(math.degrees(math.asin(throat / pitch))).tolist()


# figure 2.24
def stall_inc_negative(blade_inlet_angle, pitch, throat):
    blade_inlet_angle_range = np.arange(50, 135, 10)
    exit_angle_range = np.arange(70, 5, -10)
    i_angle_coeff = np.array([[-10, -16, -18.5, -19.5, -20, -19.5, -19, -17.5, -16],
                              [-15, -22, -23, -22, -21.5, -21, -19.5, -17.5, -16],
                              [-20.5, -28, -28.5, -27, -25, -23, -21, -18.5, -16],
                              [-26, -32.5, -33, -32, -28.5, -26, -24, -21, -18],
                              [-31, -35.5, -36, -34, -31.5, -27.5, -25, -22.5, -19.8],
                              [-35, -39, -39.3, -37, -34, -30, -27, -24, -22],
                              [-39, -42.5, -42, -40, -36, -32, -28.5, -25, -22.5]])
    f1 = [interp1d(blade_inlet_angle_range, i, fill_value='extrapolate') for i in i_angle_coeff]
    y1 = [func(blade_inlet_angle) for func in f1]
    return interp1d(exit_angle_range, y1, fill_value='extrapolate')(math.degrees(math.asin(throat / pitch))).tolist()


# figure 2.26
def stall_inc_positive(blade_inlet_angle, pitch, throat):
    blade_inlet_angle_range = np.arange(70.0, 175.0, 10.0)
    exit_angle_plus_range = np.arange(10.0, 65, 10.0)
    exit_angle_minus_range = np.arange(70.0, 15.0, -10.0)
    i_stall_plus_coeff = np.array([[2, 7, 13, 15.5, 19, 22, 23, 22.5, 18, 12, 0],
                                   [7.5, 15, 21, 25, 27, 27.1, 26.5, 21, 14, 2.5, 0],
                                   [10.5, 19, 26, 30, 30.1, 27.5, 22.5, 14, 4, 2.5, 0],
                                   [14, 21.5, 26, 26.5, 24, 20, 14, 7, 4, 2.5, 0],
                                   [16, 17.5, 17.8, 17.5, 16, 13, 8.5, 7, 4, 2.5, 0],
                                   [14, 14, 14, 13.5, 12.5, 11, 8.5, 7, 4, 2.5, 0]])
    i_stall_minus_coeff = np.array([[-18, -19.5, -20, -19.7, -18, -16.5, -16, -15, -15, -15, -15],
                                    [-23, -22.5, -22, -21, -18, -16.5, -16, -15, -15, -15, -15],
                                    [-28, -27, -25, -23, -21, -18, -16, -15, -15, -15, -15],
                                    [-32.5, -28.5, -26, -23.5, -20, -17, -16, -15, -15, -15, -15],
                                    [-36, -34, -31.5, -27, -25, -22.5, -19, -16.5, -15, -15, -15],
                                    [-39, -37.5, -34.5, -30, -27, -24, -21, -18, -16.5, -15, -15]])
    f1 = [interp1d(blade_inlet_angle_range, i, fill_value='extrapolate') for i in i_stall_plus_coeff]
    y1 = [func(blade_inlet_angle) for func in f1]
    i_stall_plus = interp1d(exit_angle_plus_range, y1, fill_value='extrapolate')(
        math.degrees(math.asin(throat / pitch))).tolist()
    f1 = [interp1d(blade_inlet_angle_range, i, fill_value='extrapolate') for i in i_stall_minus_coeff]
    y1 = [func(blade_inlet_angle) for func in f1]
    i_stall_minus = interp1d(exit_angle_minus_range, y1, fill_value='extrapolate')(
        math.degrees(math.asin(throat / pitch))).tolist()
    return i_stall_plus, i_stall_minus


# figure 2.20
def Npi(inlet_angle, outlet_angle, blade_inlet_angle, design_inc, pitch, width, throat):
    plus_inc_i_stall_pitch_width = PW_inc_positive(pitch, width, throat)
    plus_inc__i_stall_cr = CR_inc(inlet_angle, outlet_angle, pitch, width, throat)
    minus_inc_i_stall_pitch_width = PW_inc_negative(pitch, width, throat)
    if inlet_angle <= 91.0:
        plus_i_stall_basic = i_val_inlet_lower_90(blade_inlet_angle, pitch, throat)
        plus_i_stall = plus_i_stall_basic + plus_inc_i_stall_pitch_width + plus_inc__i_stall_cr

        minus_i_stall_basic = stall_inc_negative(blade_inlet_angle, pitch, throat)
        minus_i_stall = minus_i_stall_basic + minus_inc_i_stall_pitch_width
    else:
        plus_i_stall_basic, minus_i_stall_basic = stall_inc_positive(blade_inlet_angle, pitch, throat)
        plus_i_stall = plus_i_stall_basic + (
                    1 - (inlet_angle - 90) / (90 - math.degrees(math.asin(throat / pitch)))) * (
                               plus_inc_i_stall_pitch_width + plus_inc__i_stall_cr)
        minus_i_stall = minus_i_stall_basic + (1 - (inlet_angle - 90) / (
                90 - math.degrees(math.asin(throat / pitch)))) * minus_inc_i_stall_pitch_width

    Fi = FI_min(blade_inlet_angle, pitch, width)
    i_min = (plus_i_stall + Fi * minus_i_stall) / (1.0 + Fi)
    if abs(design_inc) < abs(i_min):
        I_ratio = (design_inc - i_min) / abs(minus_i_stall - i_min)
    else:
        I_ratio = (design_inc - i_min) / abs(plus_i_stall - i_min)
    I_ratio_range = np.arange(-1, 1.1, 0.5)
    Npi_coeff = np.array([2, 1.25, 1, 1.25, 2])
    return interp1d(I_ratio_range, Npi_coeff, fill_value='extrapolate')(I_ratio).tolist()


# figure 2.28
def EW_ar(width, length):
    width_length_range = np.arange(0, 5.7, 0.5)
    EW_ar_coeff = [0, 0.5, 0.94, 1.31, 1.62, 1.94, 2.28, 2.57, 2.9, 3.16, 3.45, 3.71]
    return interp1d(width_length_range, EW_ar_coeff, fill_value='extrapolate')(width / length).tolist()


# figure 2.29
def EW_bf(inet_angle, outlet_angle, pitch, width, w1, w2):
    vel_ratio_range = np.arange(0, 1.1, 0.2)
    fl_pitch_width_range = np.arange(4, 10.5, 1)
    EW_bf_factor = np.array([[3.1, 7, 10.8, 14.4, 17.9, 21.1],
                             [2.9, 6.3, 9.5, 12.7, 16, 19],
                             [2.8, 5.6, 8.7, 11.5, 14.4, 17],
                             [2.3, 5, 7.35, 9.9, 12.5, 14.8],
                             [2, 4.3, 6.4, 8.5, 10.6, 12.6],
                             [1.8, 3.5, 5.3, 7, 8.9, 10.5],
                             [1.5, 2.8, 4.15, 5.6, 7, 8.3]])
    f1 = [interp1d(vel_ratio_range, i, fill_value='extrapolate') for i in EW_bf_factor]
    y1 = [func(w1 ** 2 / w2 ** 2) for func in f1]
    fl = Fl(inet_angle, outlet_angle)
    return interp1d(fl_pitch_width_range, y1, fill_value='extrapolate')(fl * pitch / width).tolist()


def EW_loss(Npr, inet_angle, outlet_angle, pitch, width, length, w1, w2):
    ew_ar = EW_ar(width, length)
    ew_bf = EW_bf(inet_angle, outlet_angle, pitch, width, w1, w2)
    return Npr * ew_ar * ew_bf


def calculate_suction_length(inlet_angle, outlet_angle, axial_chord, throat, pitch):
    t = throat / math.cos(math.radians(outlet_angle))
    r = (axial_chord - t) / (2 * (1 - math.cos(math.radians(2 * inlet_angle)))) ** 0.5
    return 2 * math.pi * r * (inlet_angle / 360) + throat * ((1 / (throat / pitch) ** 2) - 1) ** 0.5


def calculate_profile_losses(inlet_angle, outlet_angle, pitch, width, te_radius, npr, npi):
    xp_basic = Xpb(inlet_angle, outlet_angle, pitch, width)
    te_np = Npt(te_radius, pitch, outlet_angle)
    te_dX = dXpte(te_radius, pitch)
    m_dX = dMach()
    s_dX = dBBr()
    return xp_basic * te_np * npr * npi + te_dX + m_dX + s_dX


def calculate_secondary_losses(npr, pitch, height, inlet_angle, outlet_angle, w1, w2, width):
    b_sl = EW_bf(inlet_angle, outlet_angle, pitch, width, w1, w2)
    sl_ar = EW_ar(width, height)
    return b_sl * sl_ar * npr


def calculate_losses(inlet_angle, outlet_angle, blade_inlet_angle, design_inc,
                     axial_chord, throat, pitch, te_radius, chord, height,
                     nu, ks, w1, w2):
    '''
    inlet_angle - угол входа потока
    outlet_angle - угол выхода потока
    blade_inlet_angle - скелетный вхоной угол лопатки. Если не найдешь, берем равным inlet_angle
    design_inc = 0
    axial_chord - осевая хорда, ширина профиля (B)
    throat - размер горла
    pitch - шаг
    te_radius - радиус выходной кромки
    chord - хорда
    height - высота лопатки
    nu - кинематиеская вязкость
    ks = 1e-5
    w1, w2 - скорсоти
    '''
    s_l = calculate_suction_length(inlet_angle, outlet_angle, axial_chord, throat, pitch)
    npi = Npi(inlet_angle, outlet_angle, blade_inlet_angle, design_inc, pitch, s_l, throat)
    Re = w2 * throat / nu
    npr = Npr(Re, ks, s_l)

    prof_losses = calculate_profile_losses(inlet_angle, outlet_angle, pitch, s_l, te_radius, npr, npi)
    secondary_losses = calculate_secondary_losses(npr, pitch, height, inlet_angle, outlet_angle, w1, w2, s_l)
    annulus_losses = 0.0
    return prof_losses / 100., secondary_losses / 100., annulus_losses / 100.