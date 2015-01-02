__author__ = 'Marcella Wijngaarden & Joris Schefold'
import csv
import numpy as np
import matplotlib.pyplot as plt
import os
import sys


FILE = "planets_withPandM.csv"
# FILE = "planets.csv"
FILE_NAME = "Data\\" + FILE


HEADER_ROWS = 75
BINS = 60

G = float(6.67384) * float((10**(-11)))  # Gravitational constant in  m^3 / (kg * s^2)
PI = np.pi
SOLAR_MASS = 1.9891 * float((10**30))  # Solar mass in kg
FIGURE_COUNTER = 0
def readData():
    openfile = open(FILE_NAME, 'r')
    data_reader = csv.reader(openfile)
    selected_data = []
    combination = []
    counter = 0
    weighted_occurences = {}
    for row in data_reader:
        counter += 1

        # Skip header lines
        if row[0][0] == '#' or row[0][0] == 'r':
            continue

        # Set variables to be retrieved
        distance = row[43]      # Distance (PC)
        T_star = row[51]        # T* (K)
        R_star = row[59]        # R* (R_solar)
        M_star = row[55]        # M* (M_solar)
        period = row[5]         # P (days)
        e = row[13]             # Orbital eccentricity


        if M_star != "" and period != "":
            # print "calculate distance"
            distance = getDistance(float(M_star), float(period))
        else:
            print "no distance..", M_star, period

        # Distance (PC), T* (K), R* (R_solar), M* (M_solar), P (days)
        # print distance, row[51], row[59], row[55], row[13], row[5]

        # Geometric detection correction
        probability = solarRadiusToMeters(float(R_star))/parsecToMeters(float(distance))


        # Here can be checked if crucial data is present
        if T_star == "" or R_star == "":
            print T_star, R_star, distance, "SKIP"
            continue

        # Calculate planet temperature
        if e != "":
            T_planet = getPlanetTemperature(T_star, R_star, distance, float(e))
        else:
            T_planet = getPlanetTemperature(T_star, R_star, distance)


        # Weighted occurences are corrected planet occurences
        if T_planet in weighted_occurences:
            weighted_occurences[T_planet] += 1/probability
        else:
            weighted_occurences[T_planet] = 1/probability

        combination.append([T_planet, probability])
        selected_data.append(T_planet)
        # R_star = float(solarRadiusToMeters(float(R_star_)))
        # distance = float(parsecToMeters(float(distance_)))
        # print distance, float(parsecToMeters(float(distance))), T_star, R_star, float(solarRadiusToMeters(float(R_star)))
    return selected_data, combination, weighted_occurences


def makeBarchart(weighted_occurences, N_bins):
    global FIGURE_COUNTER
    weighted_occurences_list = [0 for _ in range(N_bins)]
    temp_list = weighted_occurences.keys()
    starting_temp = 0
    temp_range = max(temp_list) - starting_temp
    temp_interval = temp_range / float(N_bins)

    for temperature in temp_list:
        weighted_occurences_list[int((temperature - starting_temp )/temp_interval)] += weighted_occurences[temperature]

    x_labels = []
    N_xticks = 6
    for i in range(N_xticks):
        x_labels.append(starting_temp + temp_range * i / N_xticks)
    x_labels.append(int(max(temp_list)))
    x_positions = np.arange(len(x_labels) - 1)*float(N_bins)/len(x_labels)
    x_positions = np.append(x_positions, N_bins)

    x_labels = []
    tick_size = 1000
    x_positions = []
    temperature = starting_temp
    while temperature < max(temp_list):
        x_labels.append(temperature)
        x_positions.append(temperature/temp_interval)
        temperature += tick_size


    plt.figure(FIGURE_COUNTER)
    FIGURE_COUNTER += 1
    plt.bar(left=range(N_bins), height=weighted_occurences_list)
    plt.tight_layout()
    plt.xticks(x_positions, x_labels)
    print x_positions, x_labels
    plt.show()
    assert False


def makeHistogram(data, combination):

    global FIGURE_COUNTER
    # Labels for plotting uncorrected R/a histogram
    TITLE_1 = "Kepler planet detection distribution"
    XLABEL_1 = "R* / a   10^5"

    # Labels for plotting planet temperature histogram
    TITLE_2 = "Kepler planet distribution"
    XLABEL_2 = "Temperature (K)"
    font = {'family' : 'normal', 'size': 18}

    plt.figure(FIGURE_COUNTER)
    FIGURE_COUNTER += 1
    plt.rc('font', **font)
    plt.title(TITLE_2)
    plt.xlabel(XLABEL_2)
    plt.ylabel("Number of planets")
    plt.hist(data, bins=30, label="Nr of planets " + str(len(data)))
    plt.legend()
    # plt.xlim(640, 740)
    # plt.ylim(0, 13)
    # print min(data), max(data)
    hist, bins = np.histogram(data, bins=BINS) #, label="Kepler planets with known distance: " + str(len(data)))  # 729
    corrections = []

    for i in range(len(bins)-1):
        probability_list = []
        for combo in combination:
            if combo[0] < bins[i+1] and combo[0] > bins[i]:
                probability_list.append(combo[1])

        if len(probability_list) != 0:
            average_bin_probability = sum(probability_list)/float(len(probability_list))
            median_bin_probability = median(probability_list)
        else:
            average_bin_probability = 1
            median_bin_probability = 1
        # print average_bin_probability
        corrections.append(median_bin_probability)

    corrected_hist = []
    corrected_hist.append(hist[0])
    for j in range(1, len(hist)):
        corrected_hist.append(float(hist[j])/float(corrections[j]))

    # hist = [float(n)/probability for n in hist]

    # Plot the resulting histogram
    center = (bins[:-1]+bins[1:])/2
    width = 0.7*(bins[1]-bins[0])

    plt.figure(FIGURE_COUNTER)
    plt.bar(center, corrected_hist, align='center', width=width)
    plt.legend()
    FIGURE_COUNTER += 1


def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]

def massFromGravitation(g, r_planet):
    """
    Takes stellar surface gravitation in Log10(cm/s^2) and planet radius in solar radii.
    Returns planet mass in kg.
    """
    g = (10**g) / 100.  # Converts units from Log10(cm/s^2) to m/s^2
    r_planet = solarRadiusToMeters(r_planet)
    return g * (r_planet**2) / G

def solarMasstoKg(mass):
    return float(mass * SOLAR_MASS)


def daysToSeconds(days):
    return float(days * (24*3600))


def metersToParsec(meters):
    return float(meters) * 3.24077929 * (10**(-17))


def solarRadiusToMeters(radius):
    return float(radius) * 6.955*(10**8)


def parsecToMeters(parsec):
    return float(parsec) * 3.08567758 * (10**16)


def getDistance(mass, period):
    """
    Use Kepler's third law to calculate the distance between planet and star given M_star and P_orbit
    """
    return metersToParsec((G * solarMasstoKg(mass) * (daysToSeconds(period)**2) / (4 * PI**2))**(1/3.))


def getPlanetTemperature(T_star_, R_star_, distance_, e=0, albedo=0):
    """
    For a given star temperature and distance calculates the temperature (K) of the planet. Based on the equality of
    L_in and L_out.
    """

    R_star = float(solarRadiusToMeters(float(R_star_)))
    distance = float(parsecToMeters(float(distance_)))
    T_star = float(T_star_)
    # return T_star * numpy.sqrt(R_star/(2*distance))
    return ((1 - albedo)*(T_star**4)*(R_star**2)/(4 * (distance**2) * (1 - e)**(1/2.)))**(1/4.)


data = readData()
makeBarchart(data[2], BINS)
makeHistogram(data[0], data[1])
plt.show()
