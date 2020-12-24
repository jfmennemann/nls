from scipy.sparse import diags


def get_D2_circulant_2nd_order(Jx, dx):
  
    D2 = diags([1, 1, -2, 1, 1], [-(Jx-1), -1, 0, 1, (Jx-1)], shape=(Jx, Jx))
    
    D2 = D2 / dx**2
    
    return D2


def get_D2_circulant_4th_order(Jx, dx):
    
    D2 = diags([16, -1, -1, 16, -30, 16, -1, -1, 16], [-(Jx-1), -(Jx-2), -2, -1, 0, 1, 2, Jx-2, Jx-1], shape=(Jx, Jx))
    
    D2 = D2 / (12 * dx**2)
    
    return D2


def get_D2_circulant_6th_order(Jx, dx):
    
    D2 = diags([270, -27, 2, 2, -27, 270, -490, 270, -27, 2, 2, -27, 270], [-(Jx-1), -(Jx-2), -(Jx-3), -3, -2, -1, 0, 1, 2, 3, Jx-3, Jx-2, Jx-1], shape=(Jx, Jx))
    
    D2 = D2 / (180 * dx**2)

    return D2


def get_D2_circulant_8th_order(Jx, dx):
    
    D2 = diags([8064, -1008, 128, -9, -9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9, -9, 128, -1008, 8064], [-(Jx-1), -(Jx-2), -(Jx-3), -(Jx-4), -4, -3, -2, -1, 0, 1, 2, 3, 4, Jx-4, Jx-3, Jx-2, Jx-1], shape=(Jx, Jx))
        
    D2 = D2 / (5040 * dx**2)
    
    return D2




