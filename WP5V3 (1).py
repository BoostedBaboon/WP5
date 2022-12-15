from math import *
import numpy as np
import itertools

# Material Changes
YieldStressList = [880, 503]
rhoList = [4430, 2810]
EModulusList = [113800, 71700]
PoissonRatioList = [0.342, 0.33]

# Output Variables
R = np.arange(0.1, 2, 0.01)
t1 = np.arange(0.001, 0.1, 0.01)
t2 = np.arange(0.001, 0.1, 0.01)
L = np.arange(0.1, 10, 0.1)

# Fixed Variables
mSpacecraft = 535
g = 9.80655
a = 6
PTank = 1.55
mSpacecraftOld = 0
VolTankMin = 0.774


# Lists
AllowedDims = []
AllowedDimsMass = []
DimensionsFin = []


# Mass Forumla
def CalulateMass(Liter, t1iter, Riter, t2iter, rho):
    V = 2 * pi * Riter * Liter * t1iter + 4 * pi * Riter ** 2 * t2iter
    Mass = rho * V
    return Mass


# Euler Buckling Formula
def EulerPress(EModulus, IEuler, AEuler, L):
    EulerCritStress = pi ** 2 * EModulus * IEuler / (AEuler * L ** 2)
    EulerSafety = EulerCritStress / 1.5

    return EulerSafety


# Shell Buckling Formula
def ShellPress(R, t1, L, PoissonRatio, PTank, EModulus):
    Lambda = sqrt(pi ** 4 * R ** 2 * t1 ** 2 / (12 * L ** 4 * (1 - PoissonRatio ** 2)))

    k = Lambda + 12 * L ** 4 * (1 - PoissonRatio ** 2) / (pi ** 4 * R ** 2 * t1 ** 2 * Lambda)
    Q = PTank / EModulus * (R / t1) ** 2

    ShellCritStress = (1.983 - 0.983 * e ** (-23.14 * Q)) * k * pi ** 2 * EModulus / (12 * (1 - PoissonRatio ** 2)) * (
                t1 / L) ** 2
    ShellSafety = ShellCritStress / 1.5

    return ShellSafety


# Min Critical stress for internal pressure
def IntPress(R, t1, t2, PTank):
    sigmacrit1 = (PTank * R) / t1
    sigmacrit2 = (PTank * R) / (2 * t2)

    sigmacritfin = min(sigmacrit1, sigmacrit2)

    sigmacritsafe = sigmacritfin * 1.5

    return sigmacritsafe

def CalcVolume(L, t1, R, t2):
    Vol = (L - 2 * R) * pi * R ** 2 + 4 / 3 * pi * R ** 3
    return Vol



for i in range(len(EModulusList)):

    YieldStress = YieldStressList[i]
    rho = rhoList[i]
    EModulus = EModulusList[i]
    PoissonRatio = PoissonRatioList[i]
    iter = [[i, j, k, m, CalulateMass(m, j, i, k, rho)] for i in R for j in t1 for k in t2 for m in
            L]  # All permutations of above variables

    if i == 0:
        print("Titanium")

    if i == 1:
        print("Aluminium")

    Loop = False
    n = 0
    while Loop is False:
        n += 1
        for i in range(len(iter)):

            Euler = False
            Shell = False
            Inter = False
            Vol = False

            Perm = iter[i]
            if Perm:
                Riter = Perm[0]
                t1iter = Perm[1]
                t2iter = Perm[2]
                Liter = Perm[3]

            else:
                break
            Iiter = pi * Riter ** 3 * t1iter
            Aiter = pi * Riter * t1iter
            StressApplied = (mSpacecraft * g) / Aiter
            # print("S App:", StressApplied)

            EulerStress = EulerPress(EModulus, Iiter, Aiter, Liter)
            if EulerStress > StressApplied:
                # print("Euler Stress,", EulerStress)
                Euler = True

            if Euler is True:

                ShellStress = ShellPress(Riter, t1iter, Liter, PoissonRatio, PTank, EModulus)
                # print("Shell Stress:", ShellStress)

                if ShellStress > StressApplied:
                    Shell = True

            if Shell is True:

                InterStress = IntPress(Riter, t1iter, t2iter, PTank)

                if InterStress < YieldStress:
                    Inter = True

            if Inter is True:
                Volume = CalcVolume(Liter, t1iter, Riter, t2iter)
                if Volume > VolTankMin:
                    Vol = True
            if Inter and Shell and Euler and Vol is True and t1iter <= Riter and t2iter <= Riter and Liter > Riter:

                AllowedDims.append(iter[i])

        for i in range(len(AllowedDims)):
            AlowedDimsIter = AllowedDims[i]

            Mass = AlowedDimsIter[4]

            AllowedDimsMass.append(Mass)

        MinMass = min(AllowedDimsMass)

        IndexMassMin = AllowedDimsMass.index(MinMass)

        DimensionsFin = AllowedDims[IndexMassMin]

        mSpacecraftOld = mSpacecraft
        mSpacecraft = 535 + DimensionsFin[4]

        print(mSpacecraftOld, mSpacecraft)
        print(abs(mSpacecraft - mSpacecraftOld))
        if abs(mSpacecraft - mSpacecraftOld) < 10:
            Loop = True
            print(mSpacecraft)
            print(abs(mSpacecraft - mSpacecraftOld))
            print()

        else:
            print(DimensionsFin)
            DimensionsFin.clear()
            AllowedDims.clear()
            AllowedDimsMass.clear()




    print("Final Dimensions are:")
    print("R = ", DimensionsFin[0])
    print("t1 = ", DimensionsFin[1])
    print("t2 = ", DimensionsFin[2])
    print("L = ", DimensionsFin[3])
    print("Mass = ", DimensionsFin[4])

    AllowedDims.clear()
    AllowedDimsMass.clear()
    DimensionsFin.clear()
    iter.clear()
