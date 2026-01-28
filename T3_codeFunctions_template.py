# The header of the corresponding code for Task 1 is relevant also here.
# It is not repeated. Please remind yourself if needed.

# Packages needed
import numpy as np
import matplotlib.pyplot as plt
# Set default font size in plots:
plt.rcParams.update({'font.size': 12})
import sys # For sys.exit()
import os # For saving plots

def calcNodePositions(nodeX, nodeY,
                      nI, nJ, pointX, pointY):
    ################################
    # DO NOT CHANGE ANYTHING HERE! #
    ################################
    # Calculates node coordinates.
    # Only changes arrays in first row of argument list!
    # Internal nodes:
    for i in range(0, nI):
        for j in range(0, nJ):
            if i > 0 and i < nI-1:
                nodeX[i,j] = 0.5*(pointX[i,0] + pointX[i-1,0])
            if j > 0 and j < nJ-1:
                nodeY[i,j] = 0.5*(pointY[0,j] + pointY[0,j-1])
    # Boundary nodes:
    nodeX[0,:]  = pointX[0,0]  # Note: corner points needed for contour plot
    nodeY[:,0]  = pointY[0,0]  # Note: corner points needed for contour plot
    nodeX[-1,:] = pointX[-1,0] # Note: corner points needed for contour plot
    nodeY[:,-1] = pointY[0,-1] # Note: corner points needed for contour plot
    
def calcDistances(dx_PE, dx_WP, dy_PN, dy_SP, dx_we, dy_sn,
                  nI, nJ, nodeX, nodeY, pointX, pointY):
    # Calculate distances in first line of argument list.
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE
    for i in range(1, nI-1):
        for j in range(1, nJ-1):
            dx_PE[i,j] = nodeX[i+1,j] - nodeX[i,j] # ADD CODE HERE x
            dx_WP[i,j] = nodeX[i,j] - nodeX[i-1,j] # ADD CODE HERE x
            dy_PN[i,j] = nodeY[i,j+1] - nodeY[i,j] # ADD CODE HERE x
            dy_SP[i,j] = nodeY[i,j] - nodeY[i,j-1] # ADD CODE HERE x
            dx_we[i,j] = pointX[i,0] - pointX[i-1,0] # ADD CODE HERE x
            dy_sn[i,j] = pointY[0,j] - pointY[0,j-1] # ADD CODE HERE x
    
def calcInterpolationFactors(fxe, fxw, fyn, fys,
                             nI, nJ, dx_PE, dx_WP, dy_PN, dy_SP, dx_we, dy_sn):
    # Calculate interpolation factors in first row of argument list.
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE
    for i in range(1, nI-1):
        for j in range(1, nJ-1):
            fxe[i,j] = 0.5 * dx_we[i,j] / dx_PE[i,j]
            fxw[i,j] = 0.5 * dx_we[i,j] / dx_WP[i,j]
            fyn[i,j] = 0.5 * dy_sn[i,j] / dy_PN[i,j]
            fys[i,j] = 0.5 * dy_sn[i,j] / dy_SP[i,j]

def initArrays(u, v, p, Fe, Fw, Fn, Fs):
    ################################
    # DO NOT CHANGE ANYTHING HERE! #
    ################################
    # Sets initial default velocity, pressure and face fluxes
    # The velocity and face flux is later kept to zero at walls 
    # (also in corner nodes, for contour plot).
    u[:,:] = 0
    v[:,:] = 0
    p[:,:] = 0
    # Sets initial default face fluxes.
    # Note that F is here supposed to include the multiplication with area
    # Keeping 'nan' where values are not needed!
    Fe[1:-1, 1:-1] = 0
    Fw[1:-1, 1:-1] = 0
    Fn[1:-1, 1:-1] = 0
    Fs[1:-1, 1:-1] = 0

def setInletVelocityAndFlux(u, v, Fe, Fw, Fn, Fs,
                            nI, nJ, rho, dx_we, dy_sn, nodeX, nodeY, grid_type, caseID):
    # Set inlet boundary values for velocity (u/v) and face fluxes (Fw/Fe/Fn/Fs).
    # The inlet velocity is 1 m/s normal to the inlet for all cases
    # Note that F is here supposed to include the multiplication with area
    # Inlet locations for imported coarse mesh:
    # Cases 1-5:   nodeX = 0, 1.864762 < nodeY < 2.0
    # Cases 6-10:  nodeX = 0, 1.946202 < nodeY < 2.0
    # Cases 11-25: nodeY = H, 1.057008 < nodeX < 1.736459
    # Inlet locations for imported fine mesh:
    # Cases 1-5:   nodeX = 0, 1.864762 < nodeY < 2.0
    # Cases 6-10:  nodeX = 0, 1.968983 < nodeY < 2.0
    # Cases 11-25: nodeY = H, 1.263541 < nodeX < 1.736459
    match grid_type:
        case 'coarse' | 'newCoarse':
            # ADD CODE HERE x
            xMin, xMax = 1.057008, 1.736459
            jB = nJ -1
            for i in range (nI):
                x = nodeX[i,jB]
                if (x > xMin) and (x < xMax):
                    u[i,jB] = 0
                    v[i, jB] = -1.0 
                    if i in range(1,nI-2):
                        Fn[i, nJ - 2] = -rho * dx_we[i, nJ - 2]
            #pass
        case 'fine':
            # ADD CODE HERE x
            xMin, xMax = 1.263541, 1.736459
            jB = nJ - 1
            el = 1e-6 #float pt tolerence

            for i in range(nI):
                x = nodeX[i,jB]
                if (x >= xMin - el) and (x <= xMax + el):
                        u[i,jB] = 0.0
                        v[i,jB] = -1.0
                        if i in range(1,nI-1):
                            Fn[i,nJ-2] = -rho * dx_we[i,nJ-2]
            #pass
        case _:
            sys.exit("Incorrect grid type!")

def initOutletFlux(Fe, Fw, Fn, Fs,
                   nI, nJ, rho, dx_we, dy_sn, nodeX, nodeY, grid_type, caseID):
    # Initialize outlet flux
    # The outlet fluxes need some non-zero value for the correction that
    # will be done later. To make it easy, we just set a flux
    # corresponding to 1 m/s at all outlet faces.
    # Note that F is here supposed to include the multiplication with area
    # Outlet locations for imported coarse mesh:
    # Cases 1-5:   nodeX = L, 0.0 < nodeY < 0.1352378
    # Cases 6-10:  nodeX = 0, 0.0 < nodeY < 0.03101719
    # Cases 11-20: nodeX = 0, 0.0 < nodeY < 0.122958
    #              nodeX = L, 0.0 < nodeY < 0.122958
    # Cases 21-25: nodeX = L, 0.0 < nodeY < 0.122958
    # Outlet locations for imported fine mesh:
    # Cases 1-5:   nodeX = L, 0.0 < nodeY < 0.1352378
    # Cases 6-10:  nodeX = 0, 0.0 < nodeY < 0.03101719
    # Cases 11-20: nodeX = 0, 0.0 < nodeY < 0.122958
    #              nodeX = L, 0.0 < nodeY < 0.122958
    # Cases 21-25: nodeX = L, 0.0 < nodeY < 0.122958
    # ADD CODE HERE
    match grid_type:
        case 'coarse' | 'fine' | 'newCoarse':
            # ADD CODE HERE x
            L = nodeX[nI-1,0]
            i = nI-2
            for j in range (1, nJ-1):
                y = nodeY[i,j]
                if (nodeX[i+1,j] == L) and (y > 0) and (y < 0.122958):
                    Fe[i,j] = rho * 1 * dy_sn[i,j]
            #pass
        case _:
            sys.exit("Incorrect grid type!")

def correctGlobalContinuity(Fe, Fw, Fn, Fs,
                            nI, nJ):
    # Make sure that outlet flux is equal to inlet flux.
    # Hint: Find inlets/outlets by convective flux into/out of the domain.
    # In this code we don't change the outlet flux later, which forces the
    # flux through the outlet boundary/ies
    # Note that F is here supposed to include the multiplication with area
    # ADD CODE HERE x
    iE = nI - 2   # east boundary face index in Fe
    iW = 1        # west boundary face index in Fw
    jN = nJ - 2   # north boundary face index in Fn
    jS = 1        # south boundary face index in Fs

    #inflow
    Fin  = 0.0
    Fin += (Fw[iW, 1:nJ-1][Fw[iW, 1:nJ-1] > 0.0]).sum() # west inflow: Fw > 0
    Fin += (Fs[1:nI-1, jS][Fs[1:nI-1, jS] > 0.0]).sum() # south inflow: Fs > 0
    Fin += (-Fe[iE, 1:nJ-1][Fe[iE, 1:nJ-1] < 0.0]).sum() # east inflow:  Fe < 0
    Fin += (-Fn[1:nI-1, jN][Fn[1:nI-1, jN] < 0.0]).sum() # north inflow: Fn < 0

    #outflow
    Fout  = 0.0
    Fout += (Fe[iE, 1:nJ-1][Fe[iE, 1:nJ-1] > 0.0]).sum() # east outflow:  Fe > 0
    Fout += (Fn[1:nI-1, jN][Fn[1:nI-1, jN] > 0.0]).sum() # north outflow: Fn > 0
    Fout += (-Fw[iW, 1:nJ-1][Fw[iW, 1:nJ-1] < 0.0]).sum() # west outflow:  Fw < 0
    Fout += (-Fs[1:nI-1, jS][Fs[1:nI-1, jS] < 0.0]).sum() # south outflow: Fs < 0

    if Fout == 0.0:
        return

    fac = Fin / Fout
    
    #scale only the outflow parts on each boundary
    mask = Fe[iE, 1:nJ-1] > 0.0
    Fe[iE, 1:nJ-1][mask] *= fac

    mask = Fn[1:nI-1, jN] > 0.0
    Fn[1:nI-1, jN][mask] *= fac

    mask = Fw[iW, 1:nJ-1] < 0.0
    Fw[iW, 1:nJ-1][mask] *= fac

    mask = Fs[1:nI-1, jS] < 0.0
    Fs[1:nI-1, jS][mask] *= fac
    #pass

    # Just for checking (if you want to):
    # globContErr = np.sum(Fe[nI-2,1:nJ-1]) - np.sum(Fw[1,1:nJ-1]) + \
    #               np.sum(Fn[1:nI-1,nJ-2]) - np.sum(Fs[1:nI-1,1])
    # print(globContErr)

def calcD(De, Dw, Dn, Ds,
          gamma, nI, nJ, dx_PE, dx_WP, dy_PN, dy_SP, dx_we, dy_sn):
    # Calculate diffusions conductances in first row of argument list.
    # Note that D is here supposed to include the multiplication with area
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE
    for i in range (1,nI-1):
        for j in range(1,nJ-1):
            De[i,j] = gamma * dy_sn[i,j] / dx_PE[i,j] # ADD CODE HERE x
            Dw[i,j] = gamma * dy_sn[i,j] / dx_WP[i,j] # ADD CODE HERE x
            Dn[i,j] = gamma * dx_we[i,j] / dy_PN[i,j] # ADD CODE HERE x
            Ds[i,j] = gamma * dx_we[i,j] / dy_SP[i,j] # ADD CODE HERE x

def calcMomEqCoeffs_FOU_CD(aE_uv, aW_uv, aN_uv, aS_uv, aP_uv,
                           nI, nJ, alphaUV, De, Dw, Dn, Ds,
                           Fe, Fw, Fn, Fs):
    # OPTIONAL!!! ONLY IF YOU ARE INTERESTED!
    # Calculate under-relaxed momentum equation coefficients, based on FOU_CD
    # (First-Order Upwind for convection and Central Differencing for diffusion),
    # using max-functions.
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # NEGLECT CONTINUITY ERROR IN CENTRAL COEFFICIENT! (Sp = 0)
    # ADD CODE HERE (IF YOU ARE INTERESTED TO TRY IT OUT)
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            aE_uv[i,j] = De[i, j] + max(-Fe[i, j], 0.0) # ADD CODE HERE x
            aW_uv[i,j] = Dw[i, j] + max( Fw[i, j], 0.0) # ADD CODE HERE x
            aN_uv[i,j] = Dn[i, j] + max(-Fn[i, j], 0.0) # ADD CODE HERE x
            aS_uv[i,j] = Ds[i, j] + max( Fs[i, j], 0.0) # ADD CODE HERE x
            aP_uv[i,j] = (aE_uv[i,j] + aW_uv[i,j] + aN_uv[i,j] + aS_uv[i,j]) / alphaUV # ADD CODE HERE x

def calcMomEqCoeffs_Hybrid(aE_uv, aW_uv, aN_uv, aS_uv, aP_uv,
                           nI, nJ, alphaUV, De, Dw, Dn, Ds,
                           Fe, Fw, Fn, Fs, fxe, fxw, fyn, fys):
    # Calculate under-relaxed momentum equation coefficients, based on Hybrid,
    # using max-functions for non-equidistant mesh.
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # NEGLECT CONTINUITY ERROR IN CENTRAL COEFFICIENT! (Sp = 0)
    # ADD CODE HERE
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            aE_uv[i,j] = max(0.0, -Fe[i, j], De[i, j] - fxe[i, j] * Fe[i, j]) # ADD CODE HERE x
            aW_uv[i,j] = max(0.0,  Fw[i, j], Dw[i, j] + fxw[i, j] * Fw[i, j]) # ADD CODE HERE x
            aN_uv[i,j] = max(0.0, -Fn[i, j], Dn[i, j] - fyn[i, j] * Fn[i, j]) # ADD CODE HERE x
            aS_uv[i,j] = max(0.0,  Fs[i, j], Ds[i, j] + fys[i, j] * Fs[i, j]) # ADD CODE HERE x
            aP_uv[i,j] = (aE_uv[i, j] + aW_uv[i, j] + aN_uv[i, j] + aS_uv[i, j]) / alphaUV # ADD CODE HERE x

def calcMomEqSu(Su_u, Su_v,
                nI, nJ, p, dx_WP, dx_PE, dy_SP, dy_PN, dx_we, dy_sn,
                alphaUV, aP_uv, u, v, fxe, fxw, fyn, fys):
    # Calculate under-relaxed momentum equation source terms,
    # based on linearly interpolated pressure face values.
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # NEGLECT CONTINUITY ERROR IN CENTRAL COEFFICIENT! (Sp = 0)
    # ADD CODE HERE
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            p_e = fxe[i, j] * p[i+1, j] + (1.0 - fxe[i, j]) * p[i, j]
            p_w = fxw[i, j] * p[i-1, j] + (1.0 - fxw[i, j]) * p[i, j]
            p_n = fyn[i, j] * p[i, j+1] + (1.0 - fyn[i, j]) * p[i, j]
            p_s = fys[i, j] * p[i, j-1] + (1.0 - fys[i, j]) * p[i, j]
            
            #pressure-gradient source terms
            Su_u[i,j] =  (p_w - p_e) * dy_sn[i, j] # ADD CODE HERE x
            Su_v[i,j] =  (p_s - p_n) * dx_we[i, j] # ADD CODE HERE x
            #under-relaxation source (aP_uv is already under-relaxed)
            Su_u[i, j] += (1.0 - alphaUV) * aP_uv[i, j] * u[i, j]
            Su_v[i, j] += (1.0 - alphaUV) * aP_uv[i, j] * v[i, j]

def solveGaussSeidel(phi,
                     nI, nJ, aE, aW, aN, aS, aP, Su, nLinSolIter):
    # Implement the Gauss-Seidel solver for general variable phi,
    # so it can be reused for all variables.
    # Do it in two directions, as in Task 2
    # Only change arrays in first row of argument list!
    # ADD CODE HERE
    for linSolIter in range(nLinSolIter):   
        for i in range(1,nI-1):
            for j in range(1,nJ-1):
                #pass # ADD CODE HERE
                #sweep west to east
                phi[i,j] = (aE[i,j]*phi[i+1,j] + aW[i,j]*phi[i-1,j] + aN[i,j]*phi[i,j+1] + aS[i,j]*phi[i,j-1] + Su[i,j]) / aP[i,j]
        for j in range(1,nJ-1):
            for i in range(1,nI-1):
                #pass # ADD CODE HERE
                #sweep south to north
                phi[i,j] = (aE[i,j]*phi[i+1,j] + aW[i,j]*phi[i-1,j] + aN[i,j]*phi[i,j+1] + aS[i,j]*phi[i,j-1] + Su[i,j]) / aP[i,j]
                
    #pass

def calcRhieChow_noCorr(Fe, Fw, Fn, Fs,
                        nI, nJ, rho, u, v,
                        dx_we, dy_sn, fxe, fxw, fyn, fys):
    # Calculate face fluxes for the pressure correction equation source term,
    # using central differencing of velocity on a non-equidistant mesh,
    # without Rhie & Chow correction.
    # Note that F is here supposed to include the multiplication with area
    # Easiest to implement, so start with this.
    # Gives checkerboarding - have a look!
    # DO NOT TOUCH BOUNDARY FLUXES, WHICH ARE SET WITH BOUNDARY CONDITIONS!
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE x
    for i in range(1, nI-1):
        for j in range(1, nJ-1):

            #east face flux of cell P: skip if east face is on boundary (i = nI-2)
            if i < nI-2:
                u_e = fxe[i, j] * u[i+1, j] + (1.0 - fxe[i, j]) * u[i, j]
                Fe[i, j] = rho * u_e * dy_sn[i, j]

            #west face flux of cell P: skip if west face is on boundary (i = 1)
            if i > 1:
                u_w = fxw[i, j] * u[i-1, j] + (1.0 - fxw[i, j]) * u[i, j]
                Fw[i, j] = rho * u_w * dy_sn[i, j]

            #north face flux of cell P: skip if north face is on boundary (j = nJ-2)
            if j < nJ-2:
                v_n = fyn[i, j] * v[i, j+1] + (1.0 - fyn[i, j]) * v[i, j]
                Fn[i, j] = rho * v_n * dx_we[i, j]

            #south face flux of cell P: skip if south face is on boundary (j = 1)
            if j > 1:
                v_s = fys[i, j] * v[i, j-1] + (1.0 - fys[i, j]) * v[i, j]
                Fs[i, j] = rho * v_s * dx_we[i, j]

    #pass

def calcRhieChow_equiCorr(Fe, Fw, Fn, Fs,
                          nI, nJ, rho, u, v,
                          dx_we, dy_sn, fxe, fxw, fyn, fys, aP_uv, p):
    # Calculate face fluxes for pressure correction equation source term,
    # using central differencing of velocity on a non-equidistant mesh,
    # and equidistant implementation of Rhie & Chow correction term.
    # Note that F is here supposed to include the multiplication with area
    # DO NOT TOUCH BOUNDARY FLUXES, WHICH ARE SET WITH BOUNDARY CONDITIONS!
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE x
    # Pressure-free velocities (same pressure interpolation as in calcMomEqSu)
    tiny = 1e-30

    for i in range(1, nI - 2):
        for j in range(1, nJ - 1):
            A = dy_sn[i, j]
            
            u_e_lin = fxe[i,j]*u[i+1,j] + (1.0-fxe[i,j])*u[i,j]
            aP_e = 0.5 * (aP_uv[i, j] + aP_uv[i+1, j])
            d_e = A / (aP_e + tiny)
            
            pW  = p[i-1, j]
            pP  = p[i,   j]
            pE  = p[i+1, j]
            pEE = p[i+2, j]

            u_e = u_e_lin - d_e * (pE - pP) + 0.25 * d_e * ((pEE + pE) - (pP + pW))

            flux = rho * u_e * A
            Fe[i, j] = flux
            Fw[i+1, j] = flux


    for i in range(1, nI - 1):
        for j in range(1, nJ - 2):
            A = dx_we[i, j]

            v_n_lin = fyn[i,j]*v[i,j+1] + (1.0-fyn[i,j])*v[i,j]
            aP_n = 0.5 * (aP_uv[i, j] + aP_uv[i, j+1])
            d_n = A / (aP_n + tiny)

            pS  = p[i, j-1]
            pP  = p[i, j]
            pN  = p[i, j+1]
            pNN = p[i, j+2]

            v_n = v_n_lin - d_n * (pN - pP) + 0.25 * d_n * ((pNN + pN) - (pP + pS))

            flux = rho * v_n * A
            Fn[i, j] = flux
            Fs[i, j+1] = flux
            
    #pass

def calcRhieChow_nonEquiCorr(Fe, Fw, Fn, Fs,
                             nI, nJ, rho, u, v,
                             dx_we, dy_sn, fxe, fxw, fyn, fys, aP_uv, p,
                             dx_WP, dx_PE, dy_SP, dy_PN):
    # OPTIONAL!!! ONLY IF YOU ARE INTERESTED!
    # Calculate face fluxes for pressure correction equation source term,
    # using central differencing of velocity on a non-equidistant mesh,
    # and non-equidistant implementation of Rhie & Chow correction term.
    # Note that F is here supposed to include the multiplication with area
    # DO NOT TOUCH BOUNDARY FLUXES, WHICH ARE SET WITH BOUNDARY CONDITIONS!
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE (IF YOU ARE INTERESTED TO TRY IT OUT)
    
    pass

def calcPpEqCoeffs(aE_pp, aW_pp, aN_pp, aS_pp, aP_pp, de, dw, dn, ds,
                   nI, nJ, rho, dx_we, dy_sn, fxe, fxw, fyn, fys, aP_uv):
    # Calculate pressure correction equation coefficients.
    # Make sure to treat boundary conditions correctly!
    # Note that de, dw, dn, ds are useful also later,
    # so make sure to set them and use them appropritely later.
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE x
    tiny = 1e-30

    for i in range(1, nI-1):
        for j in range(1, nJ-1):

            #interpolate aP_uv to faces
            aP_e = fxe[i, j] * aP_uv[i+1, j] + (1 - fxe[i, j]) * aP_uv[i, j]
            aP_w = fxw[i, j] * aP_uv[i-1, j] + (1 - fxw[i, j]) * aP_uv[i, j]
            aP_n = fyn[i, j] * aP_uv[i, j+1] + (1 - fyn[i, j]) * aP_uv[i, j]
            aP_s = fys[i, j] * aP_uv[i, j-1] + (1 - fys[i, j]) * aP_uv[i, j]

            #d-factors used for flux/velocity correction later
            de[i, j] = dy_sn[i, j] / (aP_e + tiny)
            dw[i, j] = dy_sn[i, j] / (aP_w + tiny)
            dn[i, j] = dx_we[i, j] / (aP_n + tiny)
            ds[i, j] = dx_we[i, j] / (aP_s + tiny)

            #pp-equation coefficients (already include area via dy_sn/dx_we)
            aE_pp[i, j] = rho * de[i, j] * dy_sn[i, j]
            aW_pp[i, j] = rho * dw[i, j] * dy_sn[i, j]
            aN_pp[i, j] = rho * dn[i, j] * dx_we[i, j]
            aS_pp[i, j] = rho * ds[i, j] * dx_we[i, j]

            #homogeneous Neumann at domain boundaries
            if i == nI-2:   # east boundary adjacent CV
                aE_pp[i, j] = 0
                de[i, j] = 0
            if i == 1:      # west boundary adjacent CV
                aW_pp[i, j] = 0
                dw[i, j] = 0
            if j == nJ-2:   # north boundary adjacent CV
                aN_pp[i, j] = 0
                dn[i, j] = 0
            if j == 1:      # south boundary adjacent CV
                aS_pp[i, j] = 0
                ds[i, j] = 0

            aP_pp[i, j] = aE_pp[i, j] + aW_pp[i, j] + aN_pp[i, j] + aS_pp[i, j]
    #pass

def calcPpEqSu(Su_pp,
               nI, nJ, Fe, Fw, Fn, Fs):
    # Calculate pressure correction equation source term
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE x
    for i in range(1, nI-1):
        for j in range(1, nJ-1):
            Su_pp[i, j] = (Fw[i, j] - Fe[i, j]) + (Fs[i, j] - Fn[i, j])
    #pass

def fixPp(Su_pp, aP_pp,
          pRef_i, pRef_j, aE_pp, aW_pp, aN_pp, aS_pp):
    # Fix pressure by forcing pp to zero in reference node, through source terms
    # MAKES CONVERGENCE POOR, SO BETTER TO SKIP IT FOR NOW. TRY IF YOU LIKE
    # ADD CODE HERE
    pass

def solveTDMA(phi,
              nI, nJ, aE, aW, aN, aS, aP, Su,
              nLinSolIter_phi):
    # Implement the TDMA solver for general variable phi,
    # so it can be reused for all variables.
    # Do it in two directions, as in Task 2.
    # Only change arrays in first row of argument list!
    # ADD CODE HERE x
    P = np.zeros((nI, nJ))
    Q = np.zeros((nI, nJ))

    for _ in range(nLinSolIter_phi):

        # x-lines: W -> E, for each j
        for j in range(1, nJ-1):
            i = 1
            a, b, c = aP[i,j], aE[i,j], aW[i,j]
            d = aN[i,j]*phi[i,j+1] + aS[i,j]*phi[i,j-1] + Su[i,j]
            P[i,j] = b/a
            Q[i,j] = (d + c*phi[i-1,j])/a

            for i in range(2, nI-2):
                a, b, c = aP[i,j], aE[i,j], aW[i,j]
                d = aN[i,j]*phi[i,j+1] + aS[i,j]*phi[i,j-1] + Su[i,j]
                den = a - c*P[i-1,j]
                P[i,j] = b/den
                Q[i,j] = (d + c*Q[i-1,j])/den

            i = nI-2
            a, b, c = aP[i,j], aE[i,j], aW[i,j]
            d = aN[i,j]*phi[i,j+1] + aS[i,j]*phi[i,j-1] + Su[i,j]
            den = a - c*P[i-1,j]
            P[i,j] = 0.0
            Q[i,j] = (d + c*Q[i-1,j] + b*phi[i+1,j])/den

            for i in range(nI-2, 0, -1):
                phi[i,j] = Q[i,j] if i == nI-2 else P[i,j]*phi[i+1,j] + Q[i,j]

        # y-lines: S -> N, for each i
        for i in range(1, nI-1):
            j = 1
            a, b, c = aP[i,j], aN[i,j], aS[i,j]
            d = aE[i,j]*phi[i+1,j] + aW[i,j]*phi[i-1,j] + Su[i,j]
            P[i,j] = b/a
            Q[i,j] = (d + c*phi[i,j-1])/a

            for j in range(2, nJ-2):
                a, b, c = aP[i,j], aN[i,j], aS[i,j]
                d = aE[i,j]*phi[i+1,j] + aW[i,j]*phi[i-1,j] + Su[i,j]
                den = a - c*P[i,j-1]
                P[i,j] = b/den
                Q[i,j] = (d + c*Q[i,j-1])/den

            j = nJ-2
            a, b, c = aP[i,j], aN[i,j], aS[i,j]
            d = aE[i,j]*phi[i+1,j] + aW[i,j]*phi[i-1,j] + Su[i,j]
            den = a - c*P[i,j-1]
            P[i,j] = 0.0
            Q[i,j] = (d + c*Q[i,j-1] + b*phi[i,j+1])/den

            for j in range(nJ-2, 0, -1):
                phi[i,j] = Q[i,j] if j == nJ-2 else P[i,j]*phi[i,j+1] + Q[i,j]
    #pass

def setPressureCorrectionLevel(pp,
                               nI, nJ, pRef_i, pRef_j):
    # Set pressure correction level explicitly
    # Only change arrays in first row of argument list!
    # ADD CODE HERE x
    pp_ref = pp[pRef_i, pRef_j]
    if np.isnan(pp_ref):
        return

    core = pp[1:nI-1, 1:nJ-1]
    mask = ~np.isnan(core)
    core[mask] -= pp_ref
    #pass

def correctPressureCorrectionBC(pp,
                                nI, nJ):
    # Correct pressure correction homogeneous Neumann boundary conditions
    # Only change arrays in first row of argument list!
    # ADD CODE HEREx
    # West/East boundaries (excluding corners)
    for j in range(1, nJ-1):
        pp[0, j]    = pp[1, j]
        pp[nI-1, j] = pp[nI-2, j]

    # South/North boundaries (excluding corners)
    for i in range(1, nI-1):
        pp[i, 0]    = pp[i, 1]
        pp[i, nJ-1] = pp[i, nJ-2]


    #pass

def correctPressure(p,
                    nI, nJ, alphaP, pp):
    # Correct pressure, using explicit under-relaxation
    # Only change arrays in first row of argument list!
    # ADD CODE HERE x
    for i in range(1, nI-1):
        for j in range(1, nJ-1):
            p[i, j] += alphaP * pp[i, j]

    p[0,0] = 0.5 * ((2 * p[1,0] - p[2,0]) + (2 * p[0,1] - p[0,2]))
    p[0,nJ-1] = 0.5 * ((2 * p[1,nJ-1] - p[2,nJ-1]) + (2 * p[0,nJ-2] - p[0,nJ-3]))
    p[nI-1,0] = 0.5 * ((2 * p[nI-2,0] - p[nI-3,0]) + (2 * p[nI-1,1] - p[nI-1,2]))
    p[nI-1,nJ-1] = 0.5 * ((2 * p[nI-2,nJ-1] - p[nI-3,nJ-1]) + (2 * p[nI-1,nJ-2] - p[nI-1,nJ-3]))
    #pass

def correctPressureBCandCorners(p,
                                nI, nJ, dx_PE, dx_WP, dy_PN, dy_SP):
    # Extrapolate pressure to boundaries, using constant gradient,
    # required to get correct Suu in u-mom. equation!
    # Only change arrays in first row of argument list!
    # ADD CODE HEREx
    # West/East boundaries (exclude corners)
    for j in range(1,nJ-1):
        p[0,j] = p[1,j] + (dx_WP[1,j] / dx_PE[1,j]) * (p[1,j]-p[2,j])
        p[nI-1,j] = p[nI-2,j] + (dx_PE[nI-2,j] / dx_WP[nI-2,j]) * (p[nI-2,j] - p[nI-3,j])

    # South/North boundaries (exclude corners)
    for i in range(1,nI-1):
        p[i,0] = p[i,1] + (dy_SP[i,1] / dy_PN[i,1]) * (p[i,1] - p[i,2])
        p[i,nJ-1] = p[i,nJ-2] + (dy_PN[i,nJ-2] / dy_SP[i,nJ-2]) * (p[i,nJ-2] - p[i,nJ-3])

    # Interpolate pressure to corners (kept so all do the same)
    p[0,0] = 0.5*(p[1,0]+p[0,1])
    p[nI-1,0] = 0.5*(p[nI-2,0]+p[nI-1,1])
    p[0,nJ-1] = 0.5*(p[1,nJ-1]+p[0,nJ-2])
    p[nI-1,nJ-1] = 0.5*(p[nI-2,nJ-1]+p[nI-1,nJ-2])
    
    #pass

def correctVelocity(u, v,
                    nI, nJ, fxe, fxw, fyn, fys, pp, dy_sn, dx_we, aP_uv):
    # Correct velocity components using pp solution (DO NOT TOUCH BOUNDARIES!)
    # Only change arrays in first row of argument list!
    # ADD CODE HERE x
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            ppe = fxe[i,j] * pp[i+1,j] + (1 - fxe[i,j]) * pp[i,j]
            ppw = fxw[i,j] * pp[i-1,j] + (1 - fxw[i,j]) * pp[i,j]
            ppn = fyn[i,j] * pp[i,j+1] + (1 - fyn[i,j]) * pp[i,j]
            pps = fys[i,j] * pp[i,j-1] + (1 - fys[i,j]) * pp[i,j]
            u[i,j] += dy_sn[i,j] * (ppw - ppe) / aP_uv[i,j]
            v[i,j] += dx_we[i,j] * (pps - ppn) / aP_uv[i,j]
    #pass

def correctOutletVelocity(u, v,
                          nI, nJ, rho, dx_we, dy_sn, nodeX, nodeY, grid_type, caseID):
    # Extraplate velocity at outlet (zero gradient for both u and v)
    # Outlet locations for imported coarse mesh:
    # Cases 1-5:   nodeX = L, 0.0 < nodeY < 0.1352378
    # Cases 6-10:  nodeX = 0, 0.0 < nodeY < 0.03101719
    # Cases 11-20: nodeX = 0, 0.0 < nodeY < 0.122958
    #              nodeX = L, 0.0 < nodeY < 0.122958
    # Cases 21-25: nodeX = L, 0.0 < nodeY < 0.122958
    # Outlet locations for imported fine mesh:
    # Cases 1-5:   nodeX = L, 0.0 < nodeY < 0.1352378
    # Cases 6-10:  nodeX = 0, 0.0 < nodeY < 0.03101719
    # Cases 11-20: nodeX = 0, 0.0 < nodeY < 0.122958
    #              nodeX = L, 0.0 < nodeY < 0.122958
    # Cases 21-25: nodeX = L, 0.0 < nodeY < 0.122958
    # ADD CODE HERE
    match grid_type:
        case 'coarse' | 'fine' | 'newCoarse':
            iB  = nI - 1
            iIn = nI - 2

            for j in range(1, nJ-1):
                y = nodeY[iB, j]
                if 0.0 < y < 0.122958:
                    u[iB, j] = u[iIn, j]
                    v[iB, j] = v[iIn, j]
            #pass
        case _:
            sys.exit("Incorrect grid type!")

def correctFaceFlux(Fe, Fw, Fn, Fs,
                    nI, nJ, rho, dy_sn, dx_we, de, dw, dn, ds, pp):
    # Correct face fluxes using pp solution (DO NOT TOUCH BOUNDARIES!)
    # Note that F is here supposed to include the multiplication with area
    # Only change arrays in first row of argument list!
    for i in range(1, nI-1):
        for j in range(1, nJ-1):
            Fe[i, j] += rho * dy_sn[i, j] * de[i, j] * (pp[i, j] - pp[i+1, j])
            Fw[i, j] += rho * dy_sn[i, j] * dw[i, j] * (pp[i-1, j] - pp[i, j])
            Fn[i, j] += rho * dx_we[i, j] * dn[i, j] * (pp[i, j] - pp[i, j+1])
            Fs[i, j] += rho * dx_we[i, j] * ds[i, j] * (pp[i, j-1] - pp[i, j])
    #pass

def calcNormalizedResiduals(res_u, res_v, res_c,
                            nI, nJ, iter, u, v,
                            aP_uv, aE_uv, aW_uv, aN_uv, aS_uv, Su_u, Su_v,
                            Fe, Fw, Fn, Fs):
    # Make normalization factors global so they are remembered for next call
    # to this function. Not an ideal way to do it, since they can potentially
    # be changed elsewhere, but we do it like this anyway and make sure not
    # to change them anywhere else.
    global F_uv, F_c
    # ADD CODE HERE
    
    # Compute residuals
    res_u.append(0) # U momentum residual
    res_v.append(0) # V momentum residual
    res_c.append(0) # Continuity residual/error
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            Ru = aP_uv[i, j]*u[i, j] - (aE_uv[i, j]*u[i+1, j] + aW_uv[i, j]*u[i-1, j] + aN_uv[i, j]*u[i, j+1] + aS_uv[i, j]*u[i, j-1] + Su_u[i, j])
            Rv = aP_uv[i, j]*v[i, j] - (aE_uv[i, j]*v[i+1, j] + aW_uv[i, j]*v[i-1, j] + aN_uv[i, j]*v[i, j+1] + aS_uv[i, j]*v[i, j-1] + Su_v[i, j])
            Rc = (Fw[i, j] - Fe[i, j]) + (Fs[i, j] - Fn[i, j])

            res_u[-1] += abs(Ru) # ADD CODE HEREx
            res_v[-1] += abs(Rv) # ADD CODE HEREx
            res_c[-1] += abs(Rc) # ADD CODE HEREx
    
    # Normalization with first non-normalized residual:        
    # Same normalization factor for u,v (based on largest initial)
    # Separate normalization factor for continuity residual/error
    if iter == 0:
        F_uv = max(res_u[-1],res_v[-1])
        F_c = res_c[-1]
    res_u[-1] = res_u[-1] / F_uv
    res_v[-1] = res_v[-1] / F_uv
    res_c[-1] = res_c[-1] / F_c

def createDefaultPlots(
               nI, nJ, pointX, pointY, nodeX, nodeY, pRef_i, pRef_j,
               caseID, grid_type, u, v, uTask2, vTask2, p,
               iter, res_u, res_v, res_c):
    ################################
    # DO NOT CHANGE ANYTHING HERE! #
    ################################
    # (Do not change any input arrays!)
    if not os.path.isdir('Figures'):
        os.makedirs('Figures')
    
    # Plot mesh
    plt.figure()
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Computational mesh\n(reference pressure node marked)')
    plt.axis('equal')
    plt.vlines(pointX[:,0],pointY[0,0],pointY[0,-1],colors = 'k',linestyles = 'dashed')
    plt.hlines(pointY[0,:],pointX[0,0],pointX[-1,0],colors = 'k',linestyles = 'dashed')
    plt.plot(nodeX, nodeY, 'ro')
    plt.plot(nodeX[pRef_i,pRef_j], nodeY[pRef_i,pRef_j], 'bo')
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_mesh.png')
    plt.show()
    
    # Plot velocity vectors
    plt.figure()
    plt.quiver(nodeX.T, nodeY.T, u.T, v.T)
    plt.title('Velocity vectors')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_velocityVectors.png')
    plt.show()
    
    # Plot u-velocity contour
    plt.figure()
    tempmap=plt.contourf(nodeX.T,nodeY.T,u.T,cmap='coolwarm',levels=30)
    cbar=plt.colorbar(tempmap)
    cbar.set_label('$[m/s]$')
    plt.title('U velocity $[m/s]$')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_uVelocityContour.png')
    plt.show()
    
    # Plot v-velocity contour
    plt.figure()
    tempmap=plt.contourf(nodeX.T,nodeY.T,v.T,cmap='coolwarm',levels=30)
    cbar=plt.colorbar(tempmap)
    cbar.set_label('$[m/s]$')
    plt.title('V velocity $[m/s]$')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_vVelocityContour.png')
    plt.show()
    # Plot pressure contour
    plt.figure()
    tempmap=plt.contourf(nodeX.T,nodeY.T,p.T,cmap='coolwarm',levels=30)
    cbar=plt.colorbar(tempmap)
    cbar.set_label('$[Pa]$')
    plt.plot(nodeX[pRef_i,pRef_j], nodeY[pRef_i,pRef_j], 'bo')
    plt.title('Pressure $[Pa]$\n(reference pressure node marked)')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_pressureContour.png')
    plt.show()
    
    # Plot velocity validation
    plt.figure()
    plt.plot(nodeX[:,int(nJ/2)], u[:,int(nJ/2)], 'b', label = 'U')
    plt.plot(nodeX[:,int(nJ/2)], uTask2[:,int(nJ/2)], 'bx', label = 'U ref')
    plt.plot(nodeX[:,int(nJ/2)], v[:,int(nJ/2)], 'r', label = 'V')
    plt.plot(nodeX[:,int(nJ/2)], vTask2[:,int(nJ/2)], 'rx', label = 'V ref')
    plt.title('Velocity validation (horizontal centerline)')
    plt.xlabel('x [m]')
    plt.ylabel('Velocity [m/s]')
    plt.legend()
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_velocityValidation.png')
    plt.show()
    
    # Plot residual convergence
    plt.figure()
    plt.semilogy(range(0,iter+1), res_u, 'blue', label = 'U')
    plt.semilogy(range(0,iter+1), res_v, 'red', label = 'V')
    plt.semilogy(range(0,iter+1), res_c, 'green', label = 'Continuity')
    plt.title('Residual convergence')
    plt.xlabel('Iterations')
    plt.ylabel('Residual [-]')
    plt.legend()
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_residualConvergence.png')
    plt.show()

def createAdditionalPlots():
    pass

