__author__ = 'Marcella Wijngaarden & Joris Schefold'
import csv
import numpy as np
import matplotlib.pyplot as plt
import copy
import time
import scipy

# FILE = "planets_withPandM.csv"
# FILE = "planets.csv"
FILE = "cumulative.csv"
FILE_NAME = "Data\\" + FILE

BINS = 60
G = float(6.67384) * float((10**(-11)))  # Gravitational constant in  m^3 / (kg * s^2)
PI = np.pi
SOLAR_MASS = 1.9891 * float((10**30))  # Solar mass in kg
FIGURE_COUNTER = 0
YEAR = 365.25
MISSIONLENGTH = 4*YEAR


class Planet:
    def __init__(self, row, label_row):
        # Set variables to be retrieved
        self.distance = row[label_row.index("koi_sma")]      # Distance (AU)
        self.T_star = row[label_row.index("koi_steff")]        # T* (K)
        self.R_star = row[label_row.index("koi_srad")]        # R* (R_solar)
        self.M_star = row[label_row.index("koi_smass")]        # M* (M_solar)
        self.period = float(row[label_row.index("koi_period")])         # P (days)
        self.R_planet = row[label_row.index("koi_prad")]         # Planet radius in Earth Radius

        # try:
        #     if float(self.R_planet) > 100:
        #         print self.R_planet, "<==========="
        # except:
        #     print row[label_row.index("koi_prad")]
        self.T_planet_dat = row[label_row.index("koi_teq")]
        self.gravitation = row[label_row.index("koi_slogg")]      # Stellar surface gravitation in Log10(cm/s^2)
        self.e = row[label_row.index("koi_eccen")]             # Orbital eccentricity
        self.nr_transits = row[label_row.index("koi_num_transits")]  # Number of transits
        self.confirmed = row[label_row.index("koi_disposition")]  # CONFIRMED / CANDIDATE
        # self.M_planet = row[label_row.index("pl_massj")]  # mass in jupiter masses

    def vitalDataAvailible(self):
        """
        Return True if vital data is present, False otherwise
        """
        if self.T_star == "" or self.R_star == "" or self.distance == "":
            print self.T_star, self.R_star, self.distance, "SKIP"
            return False
        return True

    def calculateDistance(self):
        if self.M_star != "" and self.period != "":
            self.distance = getDistance(float(self.M_star), float(self.period))
        else:
            assert self.R_planet != ""
            self.M_star = massFromGravitation(self.gravitation, self.R_planet)
            self.distance = getDistance(float(self.M_star), float(self.period))

    def calculateTemperature(self):
        # Use Kepler Planet temperature if available
        if self.T_planet_dat != "":
            self.T_planet = int(float(self.T_planet_dat))
        elif self.e != "":
            self.T_planet = getPlanetTemperature(self.T_star, self.R_star, self.distance, float(self.e))
        else:
            self.T_planet = getPlanetTemperature(self.T_star, self.R_star, self.distance)

        assert self.T_planet > 1

    def calculateDensity(self):
        self.density =  jupiterMassToKg(self.M_planet)/(4.0/3*np.pi*  earthRadiusToMeters(self.R_planet)**3)

    def calculateWaterContent(self):
        """
        NOG NIET AF, MOET NUMERIEK OPGELOST WORDEN
        """
        self.R_planet = (1 + .56 * x) *  kgToEarthMass(jupiterMassToKg(self.M_planet))**(0.262 * (1 - 0.138 * x))

    def setGemetricProbability(self, geometric_probability):
        self.geometric_probability = geometric_probability

    def setTimingProbability(self, timing_probability):
        self.timing_probability = timing_probability

    def calulateTotalProbability(self):
        self.probability = self.timing_probability * self.geometric_probability
        if self.confirmed == "CANDIDATE":
            self.probability *= .9

    def standAlone(self, planets, variable="T_planet", value=20):
        # return False
        for planet in planets:
            if planet == self:
                continue
            if abs(eval("self."+variable) - eval("planet."+variable)) < value:
                # print "removing planet"
                return False
        return True


def readData():
    openfile = open(FILE_NAME, 'r')
    data_reader = csv.reader(openfile)
    selected_data = []
    combination = []
    counter = 0
    weighted_occurences = {}
    planets = []

    # Skip header lines
    row = data_reader.next()
    while row[0][0] == '#':
        row = data_reader.next()
    label_row = row

    for row in data_reader:
        # if row[label_row.index("koi_disposition")] != "CONFIRMED":
        if row[label_row.index("koi_pdisposition")] != "CANDIDATE":
            continue

        counter += 1
        planet = Planet(row, label_row)
        if not planet.vitalDataAvailible():
            continue
        planet.calculateDistance()
        planet.calculateTemperature()
        # planet.calculateDensity()

        if float(planet.R_planet) > 500:
            print "huge radius ({0}), SKIP".format(planet.R_planet)
            continue

        # if planet.density > 8000:
        #     print "denser then Iron, SKIP"
        #     continue

        # Geometric detection correction
        geometric_probability = solarRadiusToMeters(float(planet.R_star))/AUToMeters(float(planet.distance))
        planet.setGemetricProbability(geometric_probability)

        timing_probability = 1
        if 2 < MISSIONLENGTH / planet.period < 3:
            time_period = MISSIONLENGTH % planet.period  # Period in which a transit must occure
            timing_probability = time_period / planet.period
        planet.setTimingProbability(timing_probability)

        planet.calulateTotalProbability()

        planets.append(planet)


    planets = removeFlukes(planets)
    for planet in planets:
        # Weighted occurences are corrected planet occurences
        if planet.T_planet in weighted_occurences:
            weighted_occurences[planet.T_planet] += 1/planet.probability
        else:
            weighted_occurences[planet.T_planet] = 1/planet.probability

        combination.append([planet.T_planet, planet.probability])
        selected_data.append(planet.T_planet)

    print "total planets", len(planets)
    return planets, combination, weighted_occurences


def removeFlukes(planets, fluke_treshold = 1.0/10**3):
    planets_ro_remove = []
    for planet in planets:
        if planet.probability < fluke_treshold:
            print "there was a really unlikelly planet with temperature: {0} and chance {1}\
             ".format(planet.T_planet, planet.probability)
            if planet.standAlone(planets):
                print "A planet with temperature {0} was remvoed".format(planet.T_planet)
                planets_ro_remove.append(planet)

    for planet in planets_ro_remove:
        planets.remove(planet)
    return planets


def makeBarchart(weighted_occurences, N_bins):
    global FIGURE_COUNTER
    temp_list = weighted_occurences.keys()

    # Set plot range
    plot_range = [0, 3000]    # [0, max(temp_list)

    # Layout plot variables
    tick_size = 500              # x-axis tick size 1000
    font = {'size': 18}
    plt.rc('font', **font)

    weighted_occurences_list = [0 for _ in range(N_bins)]
    planet_occurences_list = [0 for _ in range(N_bins)]
    starting_temp = plot_range[0]
    temp_range = plot_range[1] - plot_range[0]  #max(temp_list) - starting_temp
    temp_interval = temp_range / float(N_bins)

    HZ_count = 0
    pl_count = 0

    for temperature in temp_list:
        if temperature <= temp_range:

            weighted_occurences_list[int((temperature - starting_temp)/float(temp_interval))] += weighted_occurences[temperature]
            planet_occurences_list[int((temperature - starting_temp)/float(temp_interval))] += 1
            # pl_count += 1
            pl_count += weighted_occurences[temperature]
            if temperature >= 273 and temperature <= 373:
                HZ_count += weighted_occurences[temperature]
        else:
            continue

    x_labels = []
    x_positions = []
    temperature = starting_temp
    while temperature <= plot_range[1] and temperature >= plot_range[0]:  #max(temp_list):
        x_labels.append(temperature)
        x_positions.append(temperature / float(temp_interval))
        temperature += tick_size

    # print weighted_occurences_list
    # print planet_occurences_list
    plt.figure(FIGURE_COUNTER)
    plt.title("Corrected Kepler planet distribution")
    plt.xlabel("Temperature (K)", fontsize='17')
    plt.ylabel("Number of planets", fontsize='17')
    FIGURE_COUNTER += 1
    plt.bar(left=range(N_bins), height=weighted_occurences_list, label="Nr of planets: " + str(int(round(pl_count)))) # "Planets in HZ : " + str(int(round(HZ_count))) + " of " + str(int(round(pl_count))))
    plt.legend(loc='upper right', fontsize='17')
    plt.tight_layout()
    plt.xticks(x_positions, x_labels)


def makeHistogram(data, combination):

    global FIGURE_COUNTER
    # Labels for plotting uncorrected R/a histogram
    TITLE_1 = "Kepler planet detection distribution"
    XLABEL_1 = "R* / a   10^5"

    # Labels for plotting planet temperature histogram
    TITLE_2 = "Kepler planet distribution"
    XLABEL_2 = "Temperature (K)"
    font = {'size': 18}

    plt.figure(FIGURE_COUNTER)
    FIGURE_COUNTER += 1
    plt.rc('font', **font)
    plt.title(TITLE_2)
    plt.xlabel(XLABEL_2)
    plt.ylabel("Number of planets")
    plt.hist(data, bins=30, label="Nr of planets " + str(len(data)))
    plt.legend()
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
    # plt.legend()
    FIGURE_COUNTER += 1


def makeScatter(planets, x_axis="T_planet", y_axis="R_planet", axis = None):
    global FIGURE_COUNTER
    x_list = []
    y_list = []
    for planet in planets:
        x_list.append(float(eval("planet."+x_axis)))
        y_list.append(float(eval("planet."+y_axis)))

    plt.figure(FIGURE_COUNTER)
    FIGURE_COUNTER += 1
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    if axis != None:
        plt.axis(axis)
    plt.scatter(x_list, y_list)


def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]


def massFromGravitation(g, r_planet):
    """
    Takes stellar surface gravitation in Log10(cm/s^2) and planet radius in Jupiter radii.
    Returns planet mass in Solar Mass units.
    """
    g = (10**float(g)) / 100.  # Converts units from Log10(cm/s^2) to m/s^2
    r_planet = earthRadiusToMeters(float(r_planet))
    return (g * (r_planet**2) / G) / SOLAR_MASS


def solarMasstoKg(mass):
    return float(mass * SOLAR_MASS)


def daysToSeconds(days):
    return float(days * (24*3600))


def metersToParsec(meters):
    return float(meters) * 3.24077929 * (10**(-17))


def earthRadiusToMeters(radius):
    return float(radius) * 6371000


def jupiterRadiusToMeters(radius):
    return float(radius) * 69911000


def jupiterMassToKg(mass):
    return float(mass) * 1.89813*10**27


def solarRadiusToMeters(radius):
    return float(radius) * 6.955*(10**8)


def parsecToMeters(parsec):
    return float(parsec) * 3.08567758 * (10**16)


def metersToAU(meters):
    return float(meters) / (1.49597870700 * 10**11)


def AUToMeters(AU):
    return float(AU) * 1.49597870700 * 10**11


def kgToEarthMass(mass):
    return mass /(5.97219*10**24)


def getDistance(mass, period):
    """
    Use Kepler's third law to calculate the distance between planet and star given M_star and P_orbit
    """
    return metersToAU((G * solarMasstoKg(mass) * (daysToSeconds(period)**2) / (4 * PI**2))**(1/3.))


def getPlanetTemperature(T_star_, R_star_, distance_, e=0, albedo=.3):
    """
    For a given star temperature and distance calculates the temperature (K) of the planet. Based on the equality of
    L_in and L_out.
    """

    R_star = float(solarRadiusToMeters(float(R_star_)))
    distance = float(AUToMeters(float(distance_)))
    T_star = float(T_star_)

    return ((1 - albedo)*(T_star**4)*(R_star**2)/(4 * (distance**2) * (1 - e)**(1/2.)))**(1/4.)


data = readData()
# makeBarchart(data[2], BINS)
makeScatter(data[0])
makeScatter(data[0], axis=[-500,3000, -5, 10])
# makeHistogram(data[0], data[1])
plt.show()