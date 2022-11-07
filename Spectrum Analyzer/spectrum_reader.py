from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from datetime import datetime

counts = []
converted_plot = []
energy_range = []
peak_data = [[],[]]
k = 0
j = 0
energy_boundaries = [100, 300]
#pikes_of_interest =  [140, 145, 228, 311, 336, 364, 415, 473, 487, 537, 756, 819, 969]
pikes_of_interest =  [128]
# 145, 228, 487, 497, 537, 667, 724, 765.8, 951, 1048; 140, 145, 228, 487, 497, 537, 667, 760, 765.8
pike_det_sens = 120
spectrum_file = open('C:/Users/Lika/Desktop/SQL/Python projects/Spectrum Analyzer/17.03.2022 № 31 на 10 см нач 12-45 1000 сек.Spe', 'r')
report = open('report.txt', 'a')
irradiation_date = '03/12/2022 11:47:30'

half_life_time=20.23*24*3600
gamma_yield=0.0277
count_efficiency=0.000443
mass_attenuation_coeff=0.9782648351706037
sample_half_height=0.0
irradiation_time=119*60
thousands_of_particles=100
show_pike=80


def gauss(x, H, A, x0, sigma):
    return H + np.abs(A) * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gauss2(x, H, H1, A, A1, x0, x01, sigma, sigma1):
    return H + np.abs(A) * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + H1 + np.abs(A1) * np.exp(-(x - x01) ** 2 / (2 * sigma1 ** 2))

def gauss_fit(x, y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma], maxfev=50000)
    return popt

def gauss_fit2(x, y, x1, y1, x0, y0):
    mean = sum(x * y) / sum(y)
    mean1 = sum(x1 * y1) / sum(y1)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    sigma1 = np.sqrt(sum(y1 * (x1 - mean) ** 2) / sum(y1))
    popt, pcov = curve_fit(gauss2, x0, y0, p0=[min(y), min(y1), max(y), max(y1), mean, mean1 , sigma, sigma1], maxfev=50000)
    return popt

def moving_average(a):
    for i in range(1, len(a)-1):
        a[i] = a[i-1]+a[i]+a[i+1]
        a[i] = a[i]/3
    return a

def interpolate(energy_range, counts):
    nods = []
    nods_energy = []
    interpolation = []
    k = 0
#    counts = moving_average(counts)
    for i in range(3):  
        nods.append(counts[i])
        nods_energy.append(energy_range[i])
    for i in range(3):  
        nods.append(counts[-3+i])
        nods_energy.append(energy_range[-3+i]) 
    for i in range(3):
        k += (nods[i+3] - nods[i])/(nods_energy[i+3] - nods_energy[i])/3
    for energy in energy_range:
        interpolation.append(k*energy + (nods[0] - k*nods_energy[0]))
    return interpolation, sum(interpolation[2:-3])

def area_cutting(peak_data, filtred_sprectrum, counts, energy_range, k, sensitivity, alignment):  
    flag1 = 0
    flag2 = 0
    boundary2 = 0
    boundary1 = 0
    number = 1
    left_or_right = 0
    if k - 2 >= 0 and np.abs(peak_data[0][k-1]-peak_data[0][k-2]) <= sensitivity:
        i = 0 
        for energy in energy_range:
            if energy == peak_data[0][k-2]:
                x1 = i
                left_or_right = -1
                number+=1
            i += 1
    else:
        i = 0 
        for energy in energy_range:
            if energy == peak_data[0][k-1]:
                x1 = i
            i += 1
    if k+1 <= len(peak_data[0]) and np.abs(peak_data[0][k-1]-peak_data[0][k]) <= sensitivity:
        i = 0 
        for energy in energy_range:
            if energy == peak_data[0][k]:
                x2 = i
                number+=1
                left_or_right = 1
            i += 1
    else:
        i = 0 
        for energy in energy_range:
            if energy == peak_data[0][k-1]:
                x2 = i
            i += 1
    if min(filtred_sprectrum[x1-10:x2+10])/15 < -3:
        step = min(filtred_sprectrum[x1-10:x2+10])/15
    else:
        step = -3
    for i in range(100):
        if filtred_sprectrum[x2 + i] >= step and flag1 == 1:
            if number == 2:
                boundary2 = i+8
            elif alignment == 1:
                j=0
                for energy in energy_range:
                    if energy == peak_data[0][k]:
                        x3 = j
                    j+=1
                boundary2 = np.abs((x3-x1)//2)-2
            else:
                boundary2 = i
            flag1=2
            break
        elif filtred_sprectrum[x2 + i] < step and flag1 == 0:
            flag1 = 1
    for i in range(100):
        if filtred_sprectrum[x1 - i] >= step and flag2 == 1:
            if number == 2:
                boundary1 = i+5
            elif alignment == -1:
                j=0
                for energy in energy_range:
                    if energy == peak_data[0][k-2]:
                        x3 = j
                    j+=1
                boundary1 = np.abs((x3-x1)//2)-2
            else:
                boundary1 = i-3
            flag2=2
            break
        elif filtred_sprectrum[x1 - i] < step and flag2 == 0:
            flag2 = 1      
    # print(boundary1, boundary2)
    counts = counts[x1-boundary1:(x2+boundary2)]
    energy_range = energy_range[(x1-boundary1):(x2+boundary2)] 
    filtred_sprectrum = filtred_sprectrum[x1-boundary1:(x2+boundary2)]
    result = []
    result.append(energy_range)
    result.append(counts)
    result.append(filtred_sprectrum)
    interpolation, F = interpolate(result[0], result[1])
    result.append(interpolation)
    
    for i in range(len(result[3])):
        result[3][i] = result[1][i] - result[3][i]
    if min(result[3]) < 0:
        minimum = min(result[3])
        for i in range(len(result[3])):
            result[3][i] -= minimum
    return result, number, left_or_right, F

def cut_one_pike(peak_data, filtred_sprectrum, counts, energy_range, k):  
    boundary2 = 20
    boundary1 = 20
    i = 0 
    for energy in energy_range:
            if energy == peak_data[0][k-1]:
                x = i
            i += 1
    x = counts.index(max(counts[x-3:x+3]))
    A = (counts[x-1]-counts[x])/(energy_range[x-1]-energy_range[x])
    B = counts[x]-A*energy_range[x]
    for i in range(2, 20):
#        print(energy_range[x-i], counts[x-i], (A*energy_range[x-i] + B))
        if counts[x-i] > (A*energy_range[x-i] + B):
            boundary1 = i+1
            break
    A = (counts[x]-counts[x+1])/(energy_range[x]-energy_range[x+1])
    B = counts[x]-A*energy_range[x]
    for i in range(2, 20):
        if counts[x+i] > (A*energy_range[x+i] + B):
            boundary2 = i+1
            break
#    print(boundary2,boundary1)
    counts = counts[x-boundary1:(x+boundary2)]
    energy_range = energy_range[(x-boundary1):(x+boundary2)] 
    filtred_sprectrum = filtred_sprectrum[x-boundary1:(x+boundary2)]
    result = []
    result.append(energy_range)
    result.append(counts)
    result.append(filtred_sprectrum)
    interpolation, F = interpolate(result[0], result[1])
    result.append(interpolation)
    
    for i in range(len(result[3])):
        result[3][i] = result[1][i] - result[3][i]
    if min(result[3]) < 0:
        minimum = min(result[3])
        for i in range(len(result[3])):
            result[3][i] -= minimum
    return result

def activity_calculation(half_life_time, mass_attenuation_coeff,
                         sample_half_height, gamma_yield, count_time,
                         dead_time, cps, count_efficiency,
                         time_after_irradiation, irradiation_time
                         ):
    decay_constant=np.log(2)/half_life_time
    absorbtion_coeff = np.exp(-mass_attenuation_coeff*sample_half_height)
    decay_coeff = float(decay_constant*count_time/(1-np.exp(-decay_constant*count_time))*(np.exp(decay_constant*dead_time)-1)/(decay_constant*dead_time))
    activity = float(cps/(count_efficiency*absorbtion_coeff*gamma_yield)*decay_coeff*np.exp(decay_constant*time_after_irradiation))
    element_yield = float(activity/(1-np.exp(-decay_constant*irradiation_time))*(1-np.exp(-decay_constant*3600))*360/thousands_of_particles)
    return activity, element_yield
########### считывание спектра ### ########
for line in spectrum_file:
    if j == 1:
        MEAS_TIM = np.float64(line.split(' ')[0])
        DEAD_TIM = np.float64(line.split(' ')[1])
        j = 0
    elif j == 2:
        DATE_MEA = line
        DATE_MEA = DATE_MEA[0:-1]
        j = 0
    if line == '$DATE_MEA:\n':
        j = 2
    if line == '$MEAS_TIM:\n':
        j = 1
    if line == '$ROI:\n':
        k = 0
    if k == 1:
        if len(counts) == 16383:
            k = 0
        counts.append(float(line))
    if line == '0 16383\n':
        k = 1
energy_range = np.linspace(0, 1568.5, len(counts))
print(len(counts))
date_format_str = '%m/%d/%Y %H:%M:%S'
start = datetime.strptime(irradiation_date, date_format_str)
end =   datetime.strptime(DATE_MEA, date_format_str)
print(end)
diff = end - start
hours_after_irradiation = diff.total_seconds()/3600
print(round(diff.total_seconds()/3600, 2), 'hours', '/', round(diff.total_seconds()/3600/24, 2), 'days')
report.write('START OF REPORT\n')
report.write('irradiation date: ' + str(start) + '\n')
report.write('aquisition date: ' + str(end) + '\n')

########### фильтрация ###########
counts=moving_average(counts)
# counts=moving_average(counts)

filtred_sprectrum=[]
for i in range(len(counts)):
    converted_plot.append(counts[i])
for i in range(5, len(counts)-5):
    converted_plot[i] = -(counts[i-3]+counts[i+3]+counts[i-4]+counts[i+4])/4+(counts[i-1]+counts[i]+counts[i-2]+counts[i+2]+counts[i+1])/5
for i in range(len(converted_plot)):
    filtred_sprectrum.append(converted_plot[i])
for i in range(430):
    converted_plot[i] = 0
for i in range(2000):
    if converted_plot[i] < max(converted_plot)/pike_det_sens*2:
        converted_plot[i] = 0
for i in range(2000, 4000):
    if converted_plot[i] < max(converted_plot)/pike_det_sens:
        converted_plot[i] = 0
for i in range(4000, 5000):
    if converted_plot[i] < max(converted_plot)/pike_det_sens*0.65:
        converted_plot[i] = 0
for i in range(5000, 8000):
    if converted_plot[i] < max(converted_plot)/pike_det_sens*0.25:
        converted_plot[i] = 0
for i in range(8000, 8192):
    converted_plot[i] = 0     
########### Поиски пиков ###########

flag = 0
a = 0
b = 0
i = 0
for convert in converted_plot: 
    if convert > a:
        flag = 1
        a = convert
        b = energy_range[i]
    if convert == 0 and flag == 1:
        peak_data[0].append(b)
        peak_data[1].append(a)
        a = 0
        b = 0
        flag = 0
    i+=1

########### считывание из библиотеки ###########

k = 2
flag = 0
for energy in peak_data[0]:
    peak_data.append([])
    gamma_library = open('gamma_library.txt', 'r')
    pikes_number = 0
    cash = []
    for line in gamma_library:
        if float(line.split(';')[0]) >= (energy-0.1865*5) and float(line.split(';')[0]) <= (energy+0.1865*5):
            if (len(line.split(';')[2]) == 5 and float(line.split(';')[2][0:3]) <= 236):
                cash.append(line.split(';'))
                cash[pikes_number][4] = cash[pikes_number][4][0:-1]
                cash[pikes_number][1] = float(cash[pikes_number][1])
                cash[pikes_number][0] = float(cash[pikes_number][0])
                cash[pikes_number][3] = float(cash[pikes_number][3])
                flag = 1
                pikes_number += 1
            elif len(line.split(';')[2]) == 4 and ((line.split(';')[2][-1] == 'U' or line.split(';')[2][-1] == 'I') and float(line.split(';')[2][0:3]) <= 236):
                cash.append(line.split(';'))
                cash[pikes_number][4] = cash[pikes_number][4][0:-1]
                cash[pikes_number][1] = float(cash[pikes_number][1])
                cash[pikes_number][0] = float(cash[pikes_number][0])
                cash[pikes_number][3] = float(cash[pikes_number][3])
                flag = 1
                pikes_number += 1
            elif len(line.split(';')[2]) == 4 and line.split(';')[2][-1] != 'U' and line.split(';')[2][-1] != 'I':
                cash.append(line.split(';'))
                cash[pikes_number][4] = cash[pikes_number][4][0:-1]
                cash[pikes_number][1] = float(cash[pikes_number][1])
                cash[pikes_number][0] = float(cash[pikes_number][0])
                cash[pikes_number][3] = float(cash[pikes_number][3])
                flag = 1
                pikes_number += 1
            elif len(line.split(';')[2]) < 4:
                cash.append(line.split(';'))
                cash[pikes_number][4] = cash[pikes_number][4][0:-1]
                cash[pikes_number][1] = float(cash[pikes_number][1])
                cash[pikes_number][0] = float(cash[pikes_number][0])
                cash[pikes_number][3] = float(cash[pikes_number][3])
                flag = 1
                pikes_number += 1
    if flag == 0:
        peak_data[k].append('null')
    else:
        if pikes_number >= 2:
            criteria = []
            for i in range(len(cash)):
                if cash[i][1] > 100:
                    cash[i][1]=cash[i][1]/10
                if cash[i][-1] == 'h':
                    criteria.append(cash[i][1]/(cash[i][3])*2**(-1*hours_after_irradiation/cash[i][3]))
                elif cash[i][-1] == 'd':
                    criteria.append(cash[i][1]/(cash[i][3]*24)*2**(-1*hours_after_irradiation/24/cash[i][3]))
                elif cash[i][-1] == 'y':
                    criteria.append(cash[i][1]/(cash[i][3]*24*365)*2**(-1*hours_after_irradiation/24/365/cash[i][3]))
                elif cash[i][-1] == 'm':
                    criteria.append(-1)
                elif cash[i][-1] == 's':
                    criteria.append(-10)
            first_max = criteria.index(max(criteria))
            criteria[criteria.index(max(criteria))] = -1000
            second_max = criteria.index(max(criteria))
            peak_data[k].append([])
            peak_data[k][0] = cash[first_max]
            peak_data[k].append([])
            peak_data[k][1] = cash[second_max]
        else:
            peak_data[k].append(cash[0])
    flag = 0
#    print(peak_data[k])
#    print(energy)
#    print('\n')
    k += 1
########### удаление лишних данных ###########

for k in range(2, len(peak_data)):
    if len(peak_data[k]) == 2:
        if peak_data[k][0][2] == peak_data[k][1][2]:
            del peak_data[k][1]
        elif peak_data[k][0][0] == peak_data[k][1][0]:
            del peak_data[k][1]
    for i in range(k+1, len(peak_data)):
                if len(peak_data[i]) == 2 and i != k:
                    if peak_data[k][0][2] == peak_data[i][0][2]:
                        del peak_data[i][1]
                        if len(peak_data[k]) == 2:
                            del peak_data[k][1]
                    elif peak_data[k][0][2] == peak_data[i][1][2]:
                        del peak_data[i][0]
                        if len(peak_data[k]) == 2:
                            del peak_data[k][1]
                if len(peak_data[i]) == 2 and i != k:
                    if peak_data[k][0][2][0:3] == peak_data[i][0][2][0:3]:
                        del peak_data[i][1]
                        if len(peak_data[k]) == 2:
                            del peak_data[k][1]
#for k in range(2, len(peak_data)):
#     if len(peak_data[k]) == 2:
#         for i in range(1, len(peak_data)):
#                     if len(peak_data[i]) == 2 and i != k:
#                         if peak_data[k][1][2] == peak_data[i][0][2]:
#                             del peak_data[i][1]
#                         elif peak_data[k][1][2] == peak_data[i][1][2]:
#                             del peak_data[i][0]
#    print(peak_data[k])
                        
########### вывод графиков ###########

# counts = moving_average(counts)
for pike_of_interest in pikes_of_interest:
    array = []
    for x in peak_data[0]:  
        array.append(np.abs(x - pike_of_interest))
    num_of_interest = array.index(min(array))
    print(peak_data[num_of_interest+2])
    result, pikes_crossed, flag, F = area_cutting(peak_data, filtred_sprectrum, counts, energy_range, num_of_interest+1, 5, 0)
    result.append([])
    if pikes_crossed == 1:
        result = cut_one_pike(peak_data, filtred_sprectrum, counts, energy_range, num_of_interest+1)
        H, A, x0, sigma = gauss_fit(result[0], result[1]) 
        t2 = gauss(result[0], *gauss_fit(result[0], result[1]))
        t1 = gauss(result[0], 0, A, x0, sigma)
        constant = 4*(1+(2*(len(result[0]))+1)/6)*(max(result[1])-max(result[3]))/max(result[1])
        print('Count rate w/ bg noise = ',round(sum(t1)/MEAS_TIM, 2), 'counts/s')
        print('standart error = ', round(np.sqrt(sum(t1)*(1+constant))/MEAS_TIM, 2), 'counts/s')
        activity, element_yield = activity_calculation(half_life_time, mass_attenuation_coeff, sample_half_height, gamma_yield, MEAS_TIM, DEAD_TIM, sum(t1)/MEAS_TIM, count_efficiency, diff.total_seconds(), irradiation_time)
        print('activity = ', round(activity, 2), 'decay/s')
        print('yield = ', round(element_yield, 2), 'Bq/mkA*h')
        report.write(str(round(sum(t1)/MEAS_TIM, 2)) + '    ' + str(round(np.sqrt(sum(t1)*(1+constant))/MEAS_TIM, 2)) + '\n')
        
        fig2 = plt.figure(figsize=(14, 12))
        plt.plot(result[0], result[1], 'red')
        plt.plot(result[0], result[2], 'purple')
        plt.plot(result[0], gauss(result[0], *gauss_fit(result[0], result[1])), '--b', label='fit')
        plt.grid()
        plt.show()
       
    elif pikes_crossed == 2:
          result0, whatever, whatever1, whatever2 = area_cutting(peak_data, filtred_sprectrum, counts, energy_range, num_of_interest+1, 1, flag)
          if flag == -1:
              result1, whatever,  whatever1, F = area_cutting(peak_data, filtred_sprectrum, counts, energy_range, num_of_interest, 1, 1)
          elif flag == 1:
              result1, whatever,  whatever1, whatever2 = area_cutting(peak_data, filtred_sprectrum, counts, energy_range, num_of_interest+2, 1, -1)
          gauss_parameters = gauss_fit2(result0[0], result0[1], result1[0], result1[1], result[0], result[1])
          if gauss_parameters[2] < 0:
              gauss_parameters[2] = -1*gauss_parameters[2]
          t1 = gauss(result[0], 0, gauss_parameters[2], gauss_parameters[4], gauss_parameters[6])
          gauss_parameters = gauss_fit2(result0[0], result0[1], result1[0], result1[1], result[0], result[1])
          t2 = gauss(result[0], gauss_parameters[0], gauss_parameters[2], gauss_parameters[4], gauss_parameters[6])
          constant = 4*(1+(2*(len(result0[0]))+1)/6)*F/sum(t1)
          print('Count rate w/ bg noise = ',round(sum(t1)/MEAS_TIM, 2),'counts/s')
          print('standart error = ', round(np.sqrt(sum(t1)*(1+constant))/MEAS_TIM, 2), 'counts/s')
          activity, element_yield = activity_calculation(half_life_time, mass_attenuation_coeff, sample_half_height, gamma_yield, MEAS_TIM, DEAD_TIM, sum(t1)/MEAS_TIM, count_efficiency, diff.total_seconds(), irradiation_time)
          print('activity = ', round(activity, 2), 'decay/s')
          print('yield = ', round(element_yield, 2), 'Bq/mkA*h')
          report.write(str(round(sum(t1)/MEAS_TIM, 2)) + '    ' + str(round(np.sqrt(sum(t1)*(1+constant))/MEAS_TIM, 2)) + '\n')
          
          
          fig2 = plt.figure(figsize=(10, 10))
          plt.plot(result[0], result[1], 'green')
          plt.plot(result[0], gauss2(result[0], *gauss_fit2(result0[0], result0[1], result1[0], result1[1], result[0], result[1])), '--b', label='fit')
          plt.plot(result1[0], result1[1], 'orange')
          plt.plot(result0[0], result0[1], 'red')
          plt.plot(result[0], result[2], 'purple')
          plt.grid()
          plt.show()

mark = [0, 0]
for bound in energy_boundaries:
    a = 1
    for i in range(len(energy_range)):
        if np.abs(bound - energy_range[i]) < a:
            a = np.abs(bound - energy_range[i])
            if bound == energy_boundaries[0]:
                mark[0] = i
            else:
                mark[1] = i
counts=counts[mark[0]:mark[1]]
# counts=moving_average(counts)
filtred_sprectrum=filtred_sprectrum[mark[0]:mark[1]]
converted_plot = converted_plot[mark[0]:mark[1]]
energy_range = energy_range[mark[0]:mark[1]]

fig1 = plt.figure(figsize=(16, 9))
# plt.plot(energy_range, counts, 'black',  linewidth=1.1)
plt.bar(energy_range, counts, width=0.35)
#plt.plot(energy_range, filtred_sprectrum, 'purple')
#plt.plot(energy_range, converted_plot, 'red')
plt.xticks(np.arange(energy_boundaries[0], energy_boundaries[1], step=10))
k = 0
for energy in energy_range:
    for i in range(len(peak_data[0])):
        if energy == peak_data[0][i]:
            x = 0
            y = 0
#            if i != 1 and np.abs(peak_data[0][i]-peak_data[0][i-1]) < 10:
#                if np.abs(peak_data[0][i-1]-peak_data[0][i-2]) < 10:
#                    y = 6
#                    x = max(counts)/4
#                else:
#                    y = 4
#                    x = max(counts)/8
            if len(peak_data[i+2]) == 2 and counts[k]>max(counts)/show_pike:
                plt.annotate(peak_data[i+2][0][2] + '/' + peak_data[i+2][1][2] + '?', xy=(energy, counts[k]+max(counts)/45), xytext=(energy+y, counts[k]+max(counts)/25 + x),  arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), fontsize=12,rotation='vertical', ha='center', va='bottom', fontweight="bold")
            elif len(peak_data[i+2]) == 1 and peak_data[i+2] != ['null'] and counts[k]>max(counts)/show_pike:
                plt.annotate(peak_data[i+2][0][2], xy=(energy, counts[k]+max(counts)/45), xytext=(energy+y, counts[k]+max(counts)/10 + x),  arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), fontsize=12,rotation='vertical', ha='center', va='bottom', fontweight="bold")
    k += 1
plt.xlabel("Energy of Gamma quants, keV", fontsize=12)
plt.ylabel("counts", fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
#plt.grid()
plt.tick_params(axis='x', which='both', bottom=True,
                top=False, labelbottom=True)
plt.tick_params(axis='y', which='both', right=False,
                left=True, labelleft=True)
for pos in ['right', 'top']:
    plt.gca().spines[pos].set_visible(False)
plt.show()
report.write('END OF REPORT\n')
report.write('\n')
report.write('\n')