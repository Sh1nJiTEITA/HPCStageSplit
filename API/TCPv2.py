# Thermal Circuit Part
from pint import get_application_registry
import seuif97
#import seuif97plus
import numpy as np

def NO_DIM():  return 0
def DIM():     return 1

def NUMBER_OF_TH_PARAM(): return 12

__all__ = ['ThPoint']

#Dictionary provides ID's of  __Parameters to work with IAPWS97(seuif97)

name_by_index: dict = {
    (0,):"p",     # Pressure
    (1,):"t",      # Temperature
    (2,):"rho",    # Density
    (3,):"v",      # Specific volume
    (4,):"h",      # Specific enthalpy
    (5,):"s",      # Specific entropy
    (6,):"cp",     # Isobaric heat capacity
    (7,):"cv",     # Isochoric heat capacity
    (8,):"a",      # Sound velocity
    (9,):"k",      # Isoentropy coefficient
    (10,):"x",      # Steam effiency (dryness)
    (11,):"kv",     # Kinecic viscosity
}

name_by_index_tex: dict = {
    (0,):"{p}",     # Pressure
    (1,):"{t}",      # Temperature
    (2,):"{\\rho}",    # Density
    (3,):"{v}",      # Specific volume
    (4,):"{h}",      # Specific enthalpy
    (5,):"{s}",      # Specific entropy
    (6,):"{c}_{p",     # Isobaric heat capacity
    (7,):"{c}_{v",     # Isochoric heat capacity
    (8,):"{a}",      # Sound velocity
    (9,):"{k}",      # Isoentropy coefficient
    (10,):"{x}",      # Steam effiency (dryness)
    (11,):"{k}_{v",     # Kinecic viscosity
}

index_by_name: dict = {
    "p":    (0,),     # Pressure
    "t":    (1,),      # Temperature
    "rho":  (2,),    # Density
    "v":    (3,),      # Specific volume
    "h":    (4,),      # Specific enthalpy
    "s":    (5,),      # Specific entropy
    "cp":   (6,),     # Isobaric heat capacity
    "cv":   (7,),     # Isochoric heat capacity
    "a":    (8,),      # Sound velocity
    "k":    (9,),      # Isoentropy coefficient
    "x":    (10,),      # Steam effiency (dryness)
    "kv":   (11,),     # Kinecic viscosity
}

id_by_index: dict = {
    (0,):0,      # Pressure
    (1,):1,      # Temperature
    (2,):2,      # Density
    (3,):3,      # Specific volume
    (4,):4,      # Specific enthalpy
    (5,):5,      # Specific entropy
    (6,):8,      # Isobaric heat capacity
    (7,):9,      # Isochoric heat capacity
    (8,):10,     # Sound velocity
    (9,):11,     # Isoentropy coefficient
    (10,):15,    # Steam effiency (dryness)
    (11,):25,    # Kinecic viscosity
}
dim_by_index_or_name: dict = {
    (0,):      "MPa",           # Pressure
    "p":        "MPa",
    (1,):      "degC",          # Temperature
    "t":        "degC",
    (2,):      "kg/m^3",        # Density
    "rho":      "kg/m^3",
    (3,):      "m^3/kg",        # Specific volume
    "v":        "m^3/kg",
    (4,):      "kJ/kg",         # Specific enthalpy
    "h":        "kJ/kg",
    (5,):      "kJ/(kg*K)",     # Specific entropy
    "s":        "kJ/(kg*K)",
    (6,):      "kJ/(kg*K)",     # Isobaric heat capacity
    "cp":       "kJ/(kg*K)",
    (7,):      "kJ/(kg*K)",     # Isochoric heat capacity
    "cv":       "kJ/(kg*K)",
    (8,):     "m/s",           # Sound velocity
    "a":        "m/s",
    (9,):     "dimensionless", # Isoentropy coefficient
    "k":        "dimensionless",
    (10,):     "dimensionless", # Steam effiency (dryness)
    "x":        "dimensionless",
    (11,):     "m^2/s",         # Kinecic viscosity
    "kv":       "m^2/s",
}

#another stuff 
def isNaN(input_value):
    if isinstance(input_value, (list,np.ndarray)):
        return 0
    if input_value != input_value:
        return 1
    else:
        return 0    
    
class ThPoint:
    __un = get_application_registry()
    __Q  = __un.Quantity
    
    def __init__(self, 
                 p   = float("NaN"),
                 t   = float("NaN"),
                 #rho = float("NaN"),
                 v   = float("NaN"),
                 h   = float("NaN"),
                 s   = float("NaN"),
                 #cp  = float("NaN"),
                 #cv  = float("NaN"),
                 #a   = float("NaN"),
                 #k   = float("NaN"),
                 x   = float("NaN"),
                 #kv  = float("NaN")
                 ):    
        
        self.__Parameters = 0
        self.Update(p=p,t=t,v=v,h=h,s=s,x=x)
         
    def con_dim(self, input_value, id):
            local_str = str(input_value)
            if (local_str.find(" ") == -1):
                return input_value
            else:
                return input_value.to(dim_by_index_or_name.get((id,))).magnitude
    
    def Update(self, 
                 p   = float("NaN"),
                 t   = float("NaN"),
                 #rho = float("NaN"), 
                 v   = float("NaN"),
                 h   = float("NaN"),
                 s   = float("NaN"),
                 #cp  = float("NaN"),
                 #cv  = float("NaN"),
                 #a   = float("NaN"),
                 #k   = float("NaN"),
                 x   = float("NaN")
                 #kv  = float("NaN")
                ):
        
        self.__Parameters = np.array(
            [
                self.con_dim(p,   0),
                self.con_dim(t,   1),
                self.con_dim(float("NaN"),   2),
                self.con_dim(v,   3),
                self.con_dim(h,   4),
                self.con_dim(s,   5),
                self.con_dim(float("NaN"),   6),
                self.con_dim(float("NaN"),   7),
                self.con_dim(float("NaN"),   8),
                self.con_dim(float("NaN"),   9),
                self.con_dim(x,  10),
                self.con_dim(float("NaN"),  11),
            ]
        )
        
        nonNaNids = []
        for index, it in np.ndenumerate(self.__Parameters):
            if not(isNaN(it)):
                nonNaNids.append(index)
                
        if len(nonNaNids) != 2:
            raise Exception("Invalid input thermal paramters (number != 2)")
                
        for index, it in np.ndenumerate(self.__Parameters):
            if isNaN(it):

                self.__Parameters[index] = self.GetFunctionByID(
                    id_by_index.get(nonNaNids[0]),
                    id_by_index.get(nonNaNids[1]),
                    id_by_index.get(index)
                )(
                    self.__Parameters[nonNaNids[0]],
                    self.__Parameters[nonNaNids[1]]
                )

    def __str__(self) -> str:
        local_str_arr = []
        for index, it in np.ndenumerate(self.__Parameters):
            local_str_arr.append("{}\t = {}\n".format(name_by_index.get(index),self.Par(index[0])))
        
        return str("".join(local_str_arr))

    def dl(self):
        local_arr = []
        for index, it in np.ndenumerate(self.__Parameters):
            local_arr.append((name_by_index_tex.get(index),self.Par(index[0])))
        return local_arr
    
    def Par(self, param_index, flag = DIM()):
        param_index = int(param_index)
        if (param_index < 0):
            raise Exception("Input index is below zero")
        
        if (flag == DIM()):
            return self.__Q(self.__Parameters[param_index], 
                            (dim_by_index_or_name.get((param_index,))))
        
        elif (flag == NO_DIM()):
            return self.__Parameters[param_index]
        
        else:
            return self.__Q(self.__Parameters[param_index], 
                            (dim_by_index_or_name.get((param_index,)))).to(flag) 
        
    def p(  self, flag = DIM()): return self.Par(0,  flag)
    def t(  self, flag = DIM()): return self.Par(1,  flag)
    def rho(self, flag = DIM()): return self.Par(2,  flag)
    def v(  self, flag = DIM()): return self.Par(3,  flag)
    def h(  self, flag = DIM()): return self.Par(4,  flag)
    def s(  self, flag = DIM()): return self.Par(5,  flag)
    def cp( self, flag = DIM()): return self.Par(6,  flag)
    def cv( self, flag = DIM()): return self.Par(7,  flag)
    def a(  self, flag = DIM()): return self.Par(8,  flag)
    def k(  self, flag = DIM()): return self.Par(9,  flag)
    def x(  self, flag = DIM()): return self.Par(10, flag)
    def kv( self, flag = DIM()): return self.Par(11, flag)
       
    def GetArrNoDim(self):
        return self.__Parameters.copy()
    
    def GetNameValueDict(self, flag = DIM()):
        return {
            "p"   : self.p(flag),
            "t"   : self.t(flag),
            "rho" : self.rho(flag),
            "v"   : self.v(flag),
            "h"   : self.h(flag),
            "s"   : self.s(flag),
            "cp"  : self.cp(flag),
            "cv"  : self.cv(flag),
            "a"   : self.a(flag),
            "k"   : self.k(flag),
            "x"   : self.x(flag),
            "kv"  : self.kv(flag)
        }
    
    def GetFunctionByID(self,in_f_id, in_s_id, out_id):
        # (p, t)
        if   ((in_f_id == 0) and (in_s_id == 1)):
            if   (out_id == 3):  return seuif97.pt2v
            elif (out_id == 4):  return seuif97.pt2h
            elif (out_id == 5):  return seuif97.pt2s
            elif (out_id == 15): return seuif97.pt2x
            
            elif (out_id == 2):  return pt2rho
            elif (out_id == 8):  return pt2cp
            elif (out_id == 9):  return pt2cv
            elif (out_id == 10): return pt2a
            elif (out_id == 11): return pt2k
            elif (out_id == 25): return pt2kv
            
        # (p, h)
        elif ((in_f_id == 0) and (in_s_id == 4)):
            if   (out_id == 3):  return seuif97.ph2v
            elif (out_id == 1):  return seuif97.ph2t
            elif (out_id == 5):  return seuif97.ph2s
            elif (out_id == 15): return seuif97.ph2x
            
            elif (out_id == 2):  return ph2rho
            elif (out_id == 8):  return ph2cp
            elif (out_id == 9):  return ph2cv
            elif (out_id == 10): return ph2a
            elif (out_id == 11): return ph2k
            elif (out_id == 25): return ph2kv
        
        # (p, s)
        elif ((in_f_id == 0) and (in_s_id == 5)):
            if   (out_id == 3):  return seuif97.ps2v
            elif (out_id == 4):  return seuif97.ps2h
            elif (out_id == 1):  return seuif97.ps2t
            elif (out_id == 15): return seuif97.ps2x
            
            elif (out_id == 2):  return ps2rho
            elif (out_id == 8):  return ps2cp
            elif (out_id == 9):  return ps2cv
            elif (out_id == 10): return ps2a
            elif (out_id == 11): return ps2k
            elif (out_id == 25): return ps2kv
        
        # (p, v)
        elif ((in_f_id == 0) and (in_s_id == 3)):
            if   (out_id == 1):  return seuif97.pv2t
            elif (out_id == 4):  return seuif97.pv2h
            elif (out_id == 5):  return seuif97.pv2s
            elif (out_id == 15): return seuif97.pv2x
            
            elif (out_id == 2):  return pv2rho
            elif (out_id == 8):  return pv2cp
            elif (out_id == 9):  return pv2cv
            elif (out_id == 10): return pv2a
            elif (out_id == 11): return pv2k
            elif (out_id == 25): return pv2kv
        
        # (p, x)
        elif ((in_f_id == 0) and (in_s_id == 15)):
            if   (out_id == 3):  return seuif97.px2v
            elif (out_id == 4):  return seuif97.px2h
            elif (out_id == 5):  return seuif97.px2s
            elif (out_id == 1): return seuif97.px2t
            
            elif (out_id == 2):  return px2rho
            elif (out_id == 8):  return px2cp
            elif (out_id == 9):  return px2cv
            elif (out_id == 10): return px2a
            elif (out_id == 11): return px2k
            elif (out_id == 25): return px2kv
        
        # (t, h)
        elif ((in_f_id == 1) and (in_s_id == 4)):
            if   (out_id == 3):  return seuif97.th2v
            elif (out_id == 0):  return seuif97.th2p
            elif (out_id == 5):  return seuif97.th2s
            elif (out_id == 15): return seuif97.th2x
            
            elif (out_id == 2):  return th2rho
            elif (out_id == 8):  return th2cp
            elif (out_id == 9):  return th2cv
            elif (out_id == 10): return th2a
            elif (out_id == 11): return th2k
            elif (out_id == 25): return th2kv
        
        # (t, s)
        elif ((in_f_id == 1) and (in_s_id == 5)):
            if   (out_id == 3):  return seuif97.ts2v
            elif (out_id == 4):  return seuif97.ts2h
            elif (out_id == 0):  return seuif97.ts2p
            elif (out_id == 15): return seuif97.ts2x
            
            elif (out_id == 2):  return ts2rho
            elif (out_id == 8):  return ts2cp
            elif (out_id == 9):  return ts2cv
            elif (out_id == 10): return ts2a
            elif (out_id == 11): return ts2k
            elif (out_id == 25): return ts2kv
        
        # (t, v)
        elif ((in_f_id == 1) and (in_s_id == 3)):
            if   (out_id == 0):  return seuif97.tv2p
            elif (out_id == 4):  return seuif97.tv2h
            elif (out_id == 5):  return seuif97.tv2s
            elif (out_id == 15): return seuif97.tv2x
            
            elif (out_id == 2):  return tv2rho
            elif (out_id == 8):  return tv2cp
            elif (out_id == 9):  return tv2cv
            elif (out_id == 10): return tv2a
            elif (out_id == 11): return tv2k
            elif (out_id == 25): return tv2kv
        
        # (t, x)
        elif ((in_f_id == 1) and (in_s_id == 15)):
            if   (out_id == 3):  return seuif97.tx2v
            elif (out_id == 4):  return seuif97.tx2h
            elif (out_id == 5):  return seuif97.tx2s
            elif (out_id == 0): return seuif97.tx2p
            
            elif (out_id == 2):  return tx2rho
            elif (out_id == 8):  return tx2cp
            elif (out_id == 9):  return tx2cv
            elif (out_id == 10): return tx2a
            elif (out_id == 11): return tx2k
            elif (out_id == 25): return tx2kv
        
        # (h, s)
        elif ((in_f_id == 4) and (in_s_id == 5)):
            if   (out_id == 3):  return seuif97.hs2v
            elif (out_id == 0):  return seuif97.hs2p
            elif (out_id == 1):  return seuif97.hs2t
            elif (out_id == 15): return seuif97.hs2x
            
            elif (out_id == 2):  return hs2rho
            elif (out_id == 8):  return hs2cp
            elif (out_id == 9):  return hs2cv
            elif (out_id == 10): return hs2a
            elif (out_id == 11): return hs2k
            elif (out_id == 25): return hs2kv
        else:
            raise Exception("Invalid input thermal  __Parameters")
    
from datetime import datetime        
def MakeThArr(
    p:np.ndarray   = float("NaN"),
    t:np.ndarray   = float("NaN"),
    v:np.ndarray   = float("NaN"),
    h:np.ndarray   = float("NaN"),
    s:np.ndarray   = float("NaN"),
    x:np.ndarray   = float("NaN")
    ):
    t0 = datetime.now()
    
    input_values = [
            p,
            t,
            float("NaN"),
            v,
            h,
            s,
            float("NaN"),
            float("NaN"),
            float("NaN"),
            float("NaN"),
            x,
            float("NaN")
        ]
    
    local_input_id = []
    for it in range(0,len(input_values)):
        if not(isNaN(input_values[it])):

            local_input_id.append(it)
            
    #print(len(local_input_id))
    local = []
    for itf,its in zip(input_values[local_input_id[0]], input_values[local_input_id[1]]):    
        if ((local_input_id[0] == 0) and (local_input_id[1] == 1)):
            local.append(ThPoint(p=itf,t=its).GetArrNoDim())
            continue
        if ((local_input_id[0] == 0) and (local_input_id[1] == 3)):
            local.append(ThPoint(p=itf,v=its).GetArrNoDim())
            continue
        if ((local_input_id[0] == 0) and (local_input_id[1] == 4)):
            local.append(ThPoint(p=itf,h=its).GetArrNoDim()) 
            continue   
        if ((local_input_id[0] == 0) and (local_input_id[1] == 5)):
            local.append(ThPoint(p=itf,s=its).GetArrNoDim())
            continue    
        if ((local_input_id[0] == 0) and (local_input_id[1] == 10)):
            local.append(ThPoint(p=itf,x=its).GetArrNoDim())
            continue    
        if ((local_input_id[0] == 1) and (local_input_id[1] == 3)):
            local.append(ThPoint(t=itf,v=its).GetArrNoDim())
            continue
        if ((local_input_id[0] == 1) and (local_input_id[1] == 3)):
            local.append(ThPoint(t=itf,v=its).GetArrNoDim())
            continue
        if ((local_input_id[0] == 1) and (local_input_id[1] == 4)):
            local.append(ThPoint(t=itf,h=its).GetArrNoDim())
            continue   
        if ((local_input_id[0] == 1) and (local_input_id[1] == 5)):
            local.append(ThPoint(t=itf,s=its).GetArrNoDim())
            continue    
        if ((local_input_id[0] == 1) and (local_input_id[1] == 10)):
            local.append(ThPoint(t=itf,x=its).GetArrNoDim())
            continue
        if ((local_input_id[0] == 4) and (local_input_id[1] == 5)):
            local.append(ThPoint(h=itf,s=its).GetArrNoDim())
            continue
            
    tend = datetime.now()
    #print("time: {}".format(tend - t0))
    return np.asarray(local) 
    #print()
    #print(local)
    
    """
Property_names_by_ID: dict = {
    -0: "p",      # Pressure
    -1: "t",      # Temperature
    2: "rho",    # Density
    -3: "v",      # Specific volume
    -4: "h",      # Specific enthalpy
    -5: "s",      # Specific entropy
    8: "cp",     # Isobaric heat capacity
    9: "cv",     # Isochoric heat capacity
    10:"a",      # Sound velocity
    11:"k",      # Isoentropy coefficient
    -15:"x",      # Steam effiency (dryness)
    25:"kv",     # Kinecic viscosity
}
"""

# (p, t)
def pt2cp(p,t):
    return seuif97.pt(p,t,8)

def pt2rho(p,t):
    return seuif97.pt(p,t,2)

def pt2cv(p,t):
    return seuif97.pt(p,t,9)

def pt2a(p,t):
    return seuif97.pt(p,t,10)

def pt2k(p,t):
    return seuif97.pt(p,t,11)

def pt2kv(p,t):
    return seuif97.pt(p,t,25)

# (p, h)

def ph2cp(p,h):
    return seuif97.ph(p,h,8)

def ph2rho(p,h):
    return seuif97.ph(p,h,2)

def ph2cv(p,h):
    return seuif97.ph(p,h,9)

def ph2a(p,h):
    return seuif97.ph(p,h,10)

def ph2k(p,h):
    return seuif97.ph(p,h,11)

def ph2kv(p,h):
    return seuif97.ph(p,h,25)


# (p, s)
def ps2cp(p,s):
    return seuif97.ps(p,s,8)

def ps2rho(p,s):
    return seuif97.ps(p,s,2)

def ps2cv(p,s):
    return seuif97.ps(p,s,9)

def ps2a(p,s):
    return seuif97.ps(p,s,10)

def ps2k(p,s):
    return seuif97.ps(p,s,11)

def ps2kv(p,s):
    return seuif97.ps(p,s,25)

# (p, v)
def pv2cp(p,v):
    return seuif97.pv(p,v,8)

def pv2rho(p,v):
    return seuif97.pv(p,v,2)

def pv2cv(p,v):
    return seuif97.pv(p,v,9)

def pv2a(p,v):
    return seuif97.pv(p,v,10)

def pv2k(p,v):
    return seuif97.pv(p,v,11)

def pv2kv(p,v):
    return seuif97.pv(p,v,25)

# (p, x)
def px2cp(p,x):
    return seuif97.px(p,x,8)

def px2rho(p,x):
    return seuif97.px(p,x,2)

def px2cv(p,x):
    return seuif97.px(p,x,9)

def px2a(p,x):
    return seuif97.px(p,x,10)

def px2k(p,x):
    return seuif97.px(p,x,11)

def px2kv(p,x):
    return seuif97.px(p,x,25)

# (t, h)
def th2cp(t,h):
    return seuif97.th(t,h,8)

def th2rho(t,h):
    return seuif97.th(t,h,2)

def th2cv(t,h):
    return seuif97.th(t,h,9)

def th2a(t,h):
    return seuif97.th(t,h,10)

def th2k(t,h):
    return seuif97.th(t,h,11)

def th2kv(t,h):
    return seuif97.th(t,h,25)

#(t, s)
def ts2cp(t,s):
    return seuif97.ts(t,s,8)

def ts2rho(t,s):
    return seuif97.ts(t,s,2)

def ts2cv(t,s):
    return seuif97.ts(t,s,9)

def ts2a(t,s):
    return seuif97.ts(t,s,10)

def ts2k(t,s):
    return seuif97.ts(t,s,11)

def ts2kv(t,s):
    return seuif97.ts(t,s,25)

#(t, v)
def tv2cp(t,v):
    return seuif97.tv(t,v,8)

def tv2rho(t,v):
    return seuif97.tv(t,v,2)

def tv2cv(t,v):
    return seuif97.tv(t,v,9)

def tv2a(t,v):
    return seuif97.tv(t,v,10)

def tv2k(t,v):
    return seuif97.tv(t,v,11)

def tv2kv(t,v):
    return seuif97.tv(t,v,25)

# (t, x)
def tx2cp(t,x):
    return seuif97.tx(t,x,8)

def tx2rho(t,x):
    return seuif97.tx(t,x,2)

def tx2cv(t,x):
    return seuif97.tx(t,x,9)

def tx2a(t,x):
    return seuif97.tx(t,x,10)

def tx2k(t,x):
    return seuif97.tx(t,x,11)

def tx2kv(t,x):
    return seuif97.tx(t,x,25)

# (h, s)
def hs2cp(h,s):
    return seuif97.hs(h,s,8)

def hs2rho(h,s):
    return seuif97.hs(h,s,2)

def hs2cv(h,s):
    return seuif97.hs(h,s,9)

def hs2a(h,s):
    return seuif97.hs(h,s,10)

def hs2k(h,s):
    return seuif97.hs(h,s,11)

def hs2kv(h,s):
    return seuif97.hs(h,s,25)