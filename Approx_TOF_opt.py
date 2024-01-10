import numpy as np
import math
from math import pi
import time
import random

n=3
N=2**n
N2=4**n
epsilon = 0.005
r_dig = 3

# Define the matrices
I = np.array([[1, 0], [0, 1]])
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])
T = np.array([[1,0],[0, np.exp((1j*pi)/4)]])

pauli_1=[I,X,Y,Z]

def tensor_product_recursive(paulis, depth):
        if depth == 1:
            return paulis
        else:
            return [np.kron(p, q) for p in paulis for q in tensor_product_recursive(paulis, depth - 1)]

def tensor_product(operators):
    result = operators[0]
    for op in operators[1:]:
        result = np.kron(result, op)
    return result

def trace_product(W_prime, P_out, P_in):
    product = np.matmul(W_prime, np.matmul(P_out, np.matmul(W_prime.conj().T, P_in)))
    return np.trace(product)

def matrix_product_recursive(S, depth):
        if depth == 1:
            return S
        else:
            return [np.matmul(p, q) for p in S for q in matrix_product_recursive(S, depth - 1)]

def matrix_neg_check(A, B, N):
    for i in range(N):
        for j in range(N):
            if A[i][j] != -B[i][j]:
                return False
    return True

def matrix_eq_check(A,B,N): #Returns 1 if equal
    for i in range(N):
      if eq == 0:
        break
      for j in range(N):
        #if (A[i][j] - B[i][j] < 10**-6) or (B[i][j] - A[i][j] < 10**-6):
        if (A[i][j] - B[i][j] == 0):
          eq = 1
        else:
          eq = 0
          break
    return eq

def generate_pauli_n(n):
    """Generate the set of Pauli matrices for n qubits."""
    # Base Pauli matrices
    I = np.array([[1, 0], [0, 1]])
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    T = np.array([[1,0],[0, np.exp((1j*pi)/4)]])
    paulis = [I, X, Y, Z]

    return tensor_product_recursive(paulis, n)

def pauli_string_to_matrix(pauli_string):
    # Convert a Pauli string (e.g., 'IXY') to its matrix representation
    pauli_dict = {'I': I, 'X': X, 'Y': Y, 'Z': Z}
    matrix = np.kron(pauli_dict[pauli_string[0]], pauli_dict[pauli_string[1]])

    for p in pauli_string[2:]:
      matrix = np.kron(matrix, pauli_dict[pauli_string[2]])

    return matrix

def quartDecompose(N_func):
    N_decompose=[0]*n
    i = n-1
    m=4
    for j in range(n):
      N_decompose[i]=int(N_func%m)
      N_func=(N_func-N_decompose[i])/4
      i=i-1
    return N_decompose

def quartInt(N_arr):  #Arr is indexed as 0123...(n-1)
    sum = 0
    for j in range(n):
      sum = sum + N_arr[n-1-j]*pow(4,j)

    return sum

def ARR_MATRIX(P):
    for j in range(n):
      if j == 0:
        mat = pauli_1[P[0]]
      else:
        mat=np.kron(pauli_1[P[j]],mat)

    return mat

def ARR_EQ(A,B,size):  #Returns 0 if equal
    eq = 0
    for i in range(size):
      if A[i] != B[i]:
        eq = 1
        break

    return eq

pauli_I = [I]
i_tensored = tensor_product_recursive(pauli_I, n)[0]
pauli_n = generate_pauli_n(n)


def CLIFF_TEST(W_prime): #Returns 1 if Clifford

  i = 0
  cliff = 1
  for P in pauli_n:
    result = np.matmul(W_prime, P)
    tr_result = np.abs(np.trace(result) / N)
    if (tr_result != 0) and (i == 0):
      if tr_result == 1:
        break
      else:
        prev_amp = tr_result
        i = i+1
    if (tr_result != 0) and (i != 0):
      if tr_result != prev_amp:
        print("Not Clifford after amplitude test.")
        cliff = 0
        return 0

  for P_out in pauli_n:
    if cliff == 0:
      print("Exiting outer loop of conjugation test")
      return 0
    if np.allclose(P_out, i_tensored) == True:
      continue

    for P_in in pauli_n:
      if np.allclose(P_in, i_tensored) == True:
        continue
      val = abs(trace_product(W_prime, P_out, P_in)) / N
      if val == 1:
        break
      if (val != 1) and (val != 0):
        print("Not Clifford after conjugation test.")
        cliff = 0
        break

  if cliff == 1:
    print("Passed Clifford test")
    return 1

b_conj = 1 - 4 * epsilon**2 + 2 * epsilon**4
print("b_conj = ",b_conj)
#b_conj = round(b_conj,r_dig)
print("b_conj after rounding  = ",b_conj)

l_conj = 2 * epsilon
print("l_conj = ",l_conj)
#l_conj = round(l_conj,r_dig)
print("l_conj after rounding = ",l_conj)

def A_CONJ(W_prime, epsilon):

    p = 1

    for P_out in pauli_n:
        if np.allclose(P_out, i_tensored) == True:
            continue
        if p == 1:
            p = 0

        prod_out = np.matmul(W_prime, np.matmul(P_out, W_prime.conj().T))
        for P_in in pauli_n:
            if np.allclose(P_in, i_tensored) == True:
                continue
            prod_in = np.matmul(prod_out, P_in)
            val0 = abs(np.trace(prod_in)) / N
            val = round(val0,r_dig)
            print("val0,val, b_conj,l_conj = ",val0,val, b_conj,l_conj)
            #b_conj = 1 - 4 * epsilon**2 + 2 * epsilon**4

            if b_conj <= val <= 1:
                p += 1
                print("satisfies 1st ineq: p = ",p)
                if p > 1:
                    return "NO"
            else:
                print("Not satisfies ineq 1")

            if l_conj < val < b_conj:
                print("Satisfies ineq 2")
                return "NO"

    return "YES"

LB1 = []
UB1 = []
UB0 = []

for M in range(1,N2+1):
  lb_1 = (1 - epsilon**2) / np.sqrt(M) - np.sqrt(M * (2 * epsilon**2 - epsilon**4))
  LB1.append(lb_1)
  ub_1 = (1 / np.sqrt(M)) + np.sqrt(M * (2 * epsilon**2 - epsilon**4))
  UB1.append(ub_1)
  ub_0 = np.sqrt(M * (2 * epsilon**2 - epsilon**4))
  UB0.append(ub_0)

print("LB1 = ",LB1)
print("UB1 = ",UB1)
print("UB0 = ",UB0)

#for M in range(0,N2):
  #LB1[M] = round(LB1[M],r_dig)
  #UB1[M] = round(UB1[M],r_dig)
  #UB0[M] = round(UB0[M],r_dig)

print("LB1 after rounding = ",LB1)
print("UB1 after rounding = ",UB1)
print("UB0 after rounding = ",UB0)

def A_DECIDE(W,U_tilde, epsilon):

  W_prime = np.matmul(np.conj(W.T), U_tilde)
  S_c = []

  for P in pauli_n:
    result = np.matmul(W_prime, P)
    tr_result = np.abs(np.trace(result) / N)
    tr_result = round(tr_result,r_dig)
    S_c.append(tr_result)

  S_c.sort(reverse=True)

  for M in range(1,4**n+1):
      S_1 = S_c[:M]
      S_0 = S_c[M:]
      #lb1 = (1 - epsilon**2) / np.sqrt(M) - np.sqrt(M * (2 * epsilon**2 - epsilon**4))
      #ub1 = (1 / np.sqrt(M)) + np.sqrt(M * (2 * epsilon**2 - epsilon**4))
      #ub0 = np.sqrt(M * (2 * epsilon**2 - epsilon**4))
      amp = 1
      for x in S_1:
        if LB1[M-1] <= x <= UB1[M-1]:
          amp = 1
        else:
          amp = 0
          break
        if amp == 1:
          for x in S_0:
            if 0 <= x <= UB0[M-1]:
              amp = 1
            else:
              amp = 0
              break
        if amp == 1:
          print("Conj begin")
          result = A_CONJ(W_prime, epsilon)
          if result == "YES":
            execTime = time.time()-start_time
            #print("Fin Utilde",U_tilde)
            #print("Fin S_c",S_c)
            #print("Fin M",M)
            #print("Fin S_1",S_1)
            #print("Fin S_0",S_0)
            return "YES"

  return "NO"

#-------------------------Sets---------------------
S1 = []
S2 = []
S3 = []
S4 = []
S5 = []
S6 = []
S7 = []
S8 = []
S9 = []
S10 = []

tof_gen = [[1, 4, 16], [1, 4, 32], [1, 4, 48], [1, 8, 16], [1, 8, 32], [1, 8, 48], [1, 12, 16], [1, 12, 32], [1, 12, 48], [1, 21, 41], [1, 21, 45], [1, 25, 37], [1, 25, 45], [1, 29, 37], [1, 29, 41], [2, 4, 16], [2, 4, 32], [2, 4, 48], [2, 8, 16], [2, 8, 32], [2, 8, 48], [2, 12, 16], [2, 12, 32], [2, 12, 48], [2, 22, 42], [2, 22, 46], [2, 26, 38], [2, 26, 46], [2, 30, 38], [2, 30, 42], [3, 4, 16], [3, 4, 32], [3, 4, 48], [3, 8, 16], [3, 8, 32], [3, 8, 48], [3, 12, 16], [3, 12, 32], [3, 12, 48], [3, 23, 43], [3, 23, 47], [3, 27, 39], [3, 27, 47], [3, 31, 39], [3, 31, 43], [4, 21, 38], [4, 21, 39], [4, 22, 37], [4, 22, 39], [4, 23, 37], [4, 23, 38], [5, 10, 21], [5, 10, 37], [5, 10, 53], [5, 11, 21], [5, 11, 37], [5, 11, 53], [5, 17, 42], [5, 17, 43], [5, 26, 33], [5, 26, 43], [5, 27, 33], [5, 27, 42], [6, 9, 22], [6, 9, 38], [6, 9, 54], [6, 11, 22], [6, 11, 38], [6, 11, 54], [6, 18, 41], [6, 18, 43], [6, 25, 34], [6, 25, 43], [6, 27, 34], [6, 27, 41], [7, 9, 23], [7, 9, 39], [7, 9, 55], [7, 10, 23], [7, 10, 39], [7, 10, 55], [7, 19, 41], [7, 19, 42], [7, 25, 35], [7, 25, 42], [7, 26, 35], [7, 26, 41], [8, 25, 42], [8, 25, 43], [8, 26, 41], [8, 26, 43], [8, 27, 41], [8, 27, 42], [9, 17, 54], [9, 17, 55], [9, 22, 55], [9, 23, 54], [10, 18, 53], [10, 18, 55], [10, 21, 55], [10, 23, 53], [11, 19, 53], [11, 19, 54], [11, 21, 54], [11, 22, 53], [12, 29, 46], [12, 29, 47], [12, 30, 45], [12, 30, 47], [12, 31, 45], [12, 31, 46], [13, 26, 43], [13, 27, 42], [13, 28, 58], [13, 28, 59], [13, 44, 58], [13, 44, 59], [14, 25, 43], [14, 27, 41], [14, 28, 57], [14, 28, 59], [14, 44, 57], [14, 44, 59], [15, 25, 42], [15, 26, 41], [15, 28, 57], [15, 28, 58], [15, 44, 57], [15, 44, 58]]
print("Size = ",len(tof_gen))
tof_gen_t = []
beta_0 = 3/4
beta_1 = 1/4
for l in range(len(tof_gen)):
  P1_arr = quartDecompose(tof_gen[l][0])
  P1 = ARR_MATRIX(P1_arr)
  #print("P1_arr, P1 = ",P1_arr,P1)
  P2_arr = quartDecompose(tof_gen[l][1])
  P2 = ARR_MATRIX(P2_arr)
  #print("P2_arr, P2 = ",P2_arr,P2)
  P3_arr = quartDecompose(tof_gen[l][2])
  P3 = ARR_MATRIX(P3_arr)
  #print("P3_arr, P3 = ",P3_arr,P3)
  G_P123 = beta_0*i_tensored + beta_1*(P1 + P2 + P3 - np.matmul(P1, P2) - np.matmul(P2, P3) - np.matmul(P3, P1) + np.matmul(np.matmul(P1, P2), P3))
  tof_gen_t.append(G_P123)

size_gtof = len(tof_gen_t)
print("size_gtof = ",size_gtof)

for l1 in range(size_gtof):
  prod1 = tof_gen_t[l1]
  S1.append(prod1)
  for l2 in range(size_gtof):
    if l2 == l1:
      continue
    prod2 = np.matmul(prod1,tof_gen_t[l2])
    S2.append(prod2)
    for l3 in range(size_gtof):
      if l3 == l2:
       continue
      prod3 = np.matmul(prod2,tof_gen_t[l3])
      S3.append(prod3)
      #for l4 in range(size_gtof):
        #if l4 == l3:
         #continue
        #prod = np.matmul(prod,tof_gen_t[l4])
        #S4.append(prod)


size_S1 = len(S1)
size_S2 = len(S2)
size_S3 = len(S3)
#size_S4 = len(S4)

print("Size of S1 = ",size_S1)
print("Size of S2 = ",size_S2)
print("Size of S3 = ",size_S3)
#print("Size of S4 = ",size_S4)

#----------------UNITARY----------------
#W = np.array([[1,0],[0, np.exp((1j*pi)/4)]])  #T gate
theta_k = pow(2,2)  #k=0 : Z, k=1: S, k=2 : T
Rz = np.array([[1,0],[0, np.exp((1j*pi)/theta_k)]])  #Rz gate
W = np.kron(I, np.kron(I,Rz))
print("W, epsilon = ",W, epsilon)


cliff_result = CLIFF_TEST(W)
if cliff_result == 1:
  print("Found Clifford")
else:
  print("Not Clifford")
  
#------------------OPT-----------------
start_time = time.time()
succ = 0
i1 = 0

while i1  < size_S3:
  i1 = i1+1
  if succ == 1:
    break
  U_1 = S3[i1]
  for i2 in range(size_S3):
    U_tilde = np.matmul(U_1, S3[i2])
    result = A_DECIDE(W,U_tilde, epsilon)
    if result == "YES":
      print("YES")
      succ = succ+1
      break


execTime = time.time()-start_time
print("Implementation time = ",execTime)




