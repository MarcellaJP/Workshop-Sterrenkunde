__author__ = 'Marcella Wijngaarden & Joris Schefold'
import csv
import numpy as np
import matplotlib.pyplot as plt

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
        self.T_planet_dat = row[label_row.index("koi_teq")]  # Kelvin
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
    weighted_occurences_data = {}
    # Skip header lines
    row = data_reader.next()
    while row[0][0] == '#':
        row = data_reader.next()
    label_row = row

    koude_planeten = 0
    for row in data_reader:
        # if row[label_row.index("koi_disposition")] != "CONFIRMED":
        #     continue

        counter += 1

        # Set variables to be retrieved
        distance = row[label_row.index("koi_sma")]      # Distance (AU)
        T_star = row[label_row.index("koi_steff")]        # T* (K)
        R_star = row[label_row.index("koi_srad")]        # R* (R_solar)
        M_star = row[label_row.index("koi_smass")]        # M* (M_solar)
        period = float(row[label_row.index("koi_period")])         # P (days)
        R_planet = row[label_row.index("koi_prad")]         # Planet radius in Earth Radius
        T_planet_dat = row[label_row.index("koi_teq")]
        gravitation = row[label_row.index("koi_slogg")]      # Stellar surface gravitation in Log10(cm/s^2)
        e = row[label_row.index("koi_eccen")]             # Orbital eccentricity
        nr_transits = row[label_row.index("koi_num_transits")]  # Number of transits
        confirmed =  row[label_row.index("koi_disposition")]  # CONFIRMED / CANDIDATE

        # if nr_transits > 0:
        #     continue

        if M_star != "" and period != "":
            distance = getDistance(float(M_star), float(period))
        else:
            if R_planet != "":
                M_star = massFromGravitation(gravitation, R_planet)
                distance = getDistance(float(M_star), float(period))

        # Here can be checked if crucial data is present
        if T_star == "" or R_star == "" or distance == "":
            print T_star, R_star, distance, "SKIP"
            continue

        # Calculate planet temperature
        if e != "":
            T_planet = getPlanetTemperature(T_star, R_star, distance, float(e))
        else:
            T_planet = getPlanetTemperature(T_star, R_star, distance)
        if T_planet < 10:
            koude_planeten += 1
            print "cold planet"
            continue
        if T_planet < 188:
            print T_planet, period/365.25

        # Geometric detection correction
        geometric_probability = solarRadiusToMeters(float(R_star))/AUToMeters(float(distance))

        timing_probability = 1
        if 2 < MISSIONLENGTH / period < 3:
            time_period = MISSIONLENGTH % period  # Period in which a transit must occure
            timing_probability = time_period / period

        probability = geometric_probability * timing_probability

        # Use Kepler Planet temperature if available
        if T_planet_dat != "":
            T_planet = int(float(T_planet_dat))

        # Weighted occurences are corrected planet occurences
        if T_planet in weighted_occurences:
            weighted_occurences[T_planet] += 1/probability
        else:
            weighted_occurences[T_planet] = 1/probability

        # if T_planet_dat != "":
        #     # print T_planet_dat, " and ", T_planet
        #     if T_planet_dat in weighted_occurences_data:
        #         weighted_occurences_data[T_planet_dat] += 1/probability
        #     else:
        #         weighted_occurences_data[T_planet_dat] = 1/probability

        combination.append([T_planet, probability])
        selected_data.append(T_planet)
        # R_star = float(solarRadiusToMeters(float(R_star_)))
        # distance = float(parsecToMeters(float(distance_)))
        # print distance, float(parsecToMeters(float(distance))), T_star, R_star, float(solarRadiusToMeters(float(R_star)))

    # print weighted_occurences
    print len(combination)
    print 'count = ', counter
    print 'total datalines : ', counter
    # print weighted_occurences
    # print weighted_occurences_data
    return selected_data, combination, weighted_occurences


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


def solarRadiusToMeters(radius):
    return float(radius) * 6.955*(10**8)


def parsecToMeters(parsec):
    return float(parsec) * 3.08567758 * (10**16)


def metersToAU(meters):
    return float(meters) / (1.49597870700 * 10**11)


def AUToMeters(AU):
    return float(AU) * 1.49597870700 * 10**11

def getDistance(mass, period):
    """
    Use Kepler's third law to calculate the distance between planet and star given M_star and P_orbit
    """
    return metersToAU((G * solarMasstoKg(mass) * (daysToSeconds(period)**2) / (4 * PI**2))**(1/3.))


def getPlanetTemperature(T_star_, R_star_, distance_, e=0, albedo=0):
    """
    For a given star temperature and distance calculates the temperature (K) of the planet. Based on the equality of
    L_in and L_out.
    """

    R_star = float(solarRadiusToMeters(float(R_star_)))
    distance = float(AUToMeters(float(distance_)))
    T_star = float(T_star_)

    return ((1 - albedo)*(T_star**4)*(R_star**2)/(4 * (distance**2) * (1 - e)**(1/2.)))**(1/4.)


data = readData()
makeBarchart(data[2], BINS)
# makeHistogram(data[0], data[1])
plt.show()