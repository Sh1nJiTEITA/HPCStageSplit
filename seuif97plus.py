import seuif97

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