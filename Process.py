__author__ = 'Marcella Wijngaarden & Joris Schefold'
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
import copy
import time
import scipy
from scipy import stats
from scipy.optimize import curve_fit
import csv

FILE = "cumulative.csv"
FILE_NAME = "Data\\" + FILE

FILE_2 = "phl_hec_all_kepler.csv"
FILE_NAME_2 = "Data\\" + FILE_2

FILE_3 = "data_phases_GH.csv"
FILE_NAME_3 = "Data\\" + FILE_3

BINS = 60
G = float(6.67384) * float((10**(-11)))  # Gravitational constant in  m^3 / (kg * s^2)
SOLAR_MASS = 1.9891 * float((10**30))  # Solar mass in kg
EARTH_RADIUS = 6371000 # in meters
FIGURE_COUNTER = 0
YEAR = 365.25
MISSIONLENGTH = 4*YEAR
GASS_CONSTANT = 8.3144621
GREENHOUSE = False
GREENHOUSE = True
EPSILON = .78

class Planet:
    def __init__(self, row, label_row, derived_row = None, derived_labels = None):
        # Set variables to be retrieved
        self.phase = ""
        self.distance = row[label_row.index("koi_sma")]      # Distance (AU)
        self.T_star = row[label_row.index("koi_steff")]        # T* (K)
        self.R_star = row[label_row.index("koi_srad")]        # R* (R_solar)
        self.M_star = row[label_row.index("koi_smass")]        # M* (M_solar)
        self.period = float(row[label_row.index("koi_period")])         # P (days)
        self.R_planet = row[label_row.index("koi_prad")]         # Planet radius in Earth Radius
        self.T_planet_dat = row[label_row.index("koi_teq")]
        self.gravitation = row[label_row.index("koi_slogg")]      # Stellar surface gravitation in Log10(cm/s^2)
        self.e = row[label_row.index("koi_eccen")]             # Orbital eccentricity
        self.nr_transits = row[label_row.index("koi_num_transits")]  # Number of transits
        self.confirmed = row[label_row.index("koi_disposition")]  # CONFIRMED / CANDIDATE
        self.composition = derived_row[derived_labels.index("P. Composition Class")]
        self.atmposphere = derived_row[derived_labels.index("P. Atmosphere Class")]
        self.M_planet = derived_row[derived_labels.index("P. Mass (EU)")]  # mass in earth masses
        self.S_ID = derived_row[derived_labels.index("S. Name")]
        self.ID = row[label_row.index("kepoi_name")]  # KOI ID
        self.S_type = derived_row[derived_labels.index("S. Type")] #G, F, K
        self.N_planets = derived_row[derived_labels.index("S. No. Planets")] #G, F, K
        self.N_planets_HZ = derived_row[derived_labels.index("S. No. Planets HZ")] #G, F, K

        # self.S_type = derived_row[derived_labels.index("S. Type")]
        # self.S_type = derived_row[derived_labels.index("S. Type")]
        # self.confirmed = row[label_row.index("koi_disposition")]
        # self.confirmed = row[label_row.index("koi_disposition")]


    def calculatePressure(self):
        if self.M_planet != "":
            self.P_planet = 101325 * float(self.M_planet)**2/(float(self.R_planet)**3)

    def vitalDataAvailible(self):
        """
        Return True if vital data is present, False otherwise
        """
        if self.T_star == "" or self.R_star == "" or self.distance == "":   # or self.M_planet == "":
            print self.T_star, self.R_star, self.distance, "SKIP"
            return False
        return True

    def calculateDistance(self):
        if self.M_star != "" and self.period != "":
            self.distance = getDistance(float(self.M_star), float(self.period))
        else:
            assert self.R_planet != ""
            self.M_star = massFromGravitation(self.gravitation, self.R_star)
            self.distance = getDistance(float(self.M_star), float(self.period))

    def calculateTemperature(self):
        # Use Kepler Planet temperature if available
        if self.T_planet_dat != "":
            self.T_planet = int(float(self.T_planet_dat))
        elif self.e != "":
            self.T_planet = getPlanetTemperature(self.T_star, self.R_star, self.distance, float(self.e))
        else:
            self.T_planet = getPlanetTemperature(self.T_star, self.R_star, self.distance)

        if GREENHOUSE:
            self.T_planet *= (1/(1-EPSILON/2))**.25
        assert self.T_planet > 10

    def calculateDensity(self):
        self.density =  EarthMassToKg(self.M_planet)/(4.0/3*np.pi * earthRadiusToMeters(self.R_planet)**3)

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
        self.probability = float(self.timing_probability * self.geometric_probability)
        if self.confirmed == "CANDIDATE":
            self.probability *= .9

    def standAlone(self, planets, variable="T_planet", value=20):
        for planet in planets:
            if planet == self:
                continue
            if abs(eval("self."+variable) - eval("planet."+variable)) < value:
                # print "removing planet"
                return False
        return True

    def setWaterPhase(self, phase):
        self.phase = phase

def readDerivedData():
    openfile = open(FILE_NAME_2, 'r')
    data_reader = csv.reader(openfile)
    header_row = data_reader.next()
    data_dict = {}
    for row in data_reader:
        number = row[2]
        name = "K" + (8-len(number))*"0"+number
        data_dict[name] = row

    return data_dict, header_row


def readPhaseData():
    openfile = open(FILE_NAME_3, 'r')
    data_reader = csv.reader(openfile)
    data_dict = {}
    for row in data_reader:
        name = row[0]
        phase = row[3]
        data_dict[name] = phase

    return data_dict


def readData(derived_data = None, derived_label_row = None, phase_dict = None, weighted_occurences_variable = "T_planet", constraints = {}):
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

        if derived_data != None:
            planet = Planet(row, label_row, derived_data[row[2]], derived_label_row)
        else:
            planet = Planet(row, label_row)
        if not planet.vitalDataAvailible():
            continue

        planet.calculateDistance()
        planet.calculateTemperature()
        if phase_dict != None:
            try:
                planet.setWaterPhase(phase_dict[planet.ID])
                print planet.phase
            except KeyError:
                pass
        if planet.M_planet != "":
            planet.calculateDensity()
            planet.calculatePressure()

        if not sattisfiesConstraints(planet, constraints):
            continue

        if float(planet.R_planet) > 15:
            print "huge radius ({0}), SKIP".format(planet.R_planet)
            continue


        # Geometric detection correction
        geometric_probability = solarRadiusToMeters(float(planet.R_star))/AUToMeters(float(planet.distance))
        planet.setGemetricProbability(geometric_probability)

        timing_probability = 1
        if MISSIONLENGTH / planet.period < 3:
            time_period = MISSIONLENGTH % planet.period  # Period in which a transit must occur
            timing_probability = time_period / planet.period
        planet.setTimingProbability(timing_probability)

        planet.calulateTotalProbability()
        planets.append(planet)

    planets = removeFlukes(planets)
    for planet in planets:
        # Weighted occurences are corrected planet occurences
        value = float(eval("planet." + weighted_occurences_variable))
        if value in weighted_occurences:
            weighted_occurences[value] += 1/planet.probability
        else:
            weighted_occurences[value] = 1/planet.probability

        combination.append([value, planet.probability])
        selected_data.append(value)

    print "total planets", len(planets)
    print "range {2} is between {0} and {1}".format(min(selected_data), max(selected_data), weighted_occurences_variable)
    return planets, combination, weighted_occurences


def removeFlukes(planets, fluke_treshold = 1.0/10**4):
    planets_ro_remove = []
    for planet in planets:
        if planet.probability < fluke_treshold:
            print "Low probability planet with T: {0} and P: {1}\
             ".format(planet.T_planet, planet.probability)
            if planet.standAlone(planets):
                print "Removed planet with temperature {0}".format(planet.T_planet)
                planets_ro_remove.append(planet)

    for planet in planets_ro_remove:
        planets.remove(planet)
    return planets


def sattisfiesConstraints(planet, constraints):
    for constraint in constraints:
        if len(constraints[constraint]) == 1:
            if constraints[constraint] != eval("planet." + constraint):
                return False
        else:
            value = eval("planet."+constraint)
            if value == "":
                return False
            else:
                value = float(value)
                if not float(constraints[constraint][0]) <= value <= float(constraints[constraint][1]):
                    return False


    return True


def getConstrainedPlanets(planets, constraints):
    """
    constraints in form of {"T_planet":(273, 373)}
    """

    total_planets = []
    for planet in planets:
        if sattisfiesConstraints(planet, constraints):
            total_planets.append(planet)

    return total_planets


def calculatePercentage(planets, constraints):
    """
    constraints in form of {"T_planet":(273, 373)}
    """

    def sattisfiesConstraints(planet, constraints):
        for constraint in constraints:
            value = float(eval("planet."+constraint))
            if not float(constraints[constraint][0]) <= value <= float(constraints[constraint][1]):
                return False
        return True

    total_planets = 0
    constrained_planets = 0
    for planet in planets:
        total_planets += 1.0/planet.probability
        if sattisfiesConstraints(planet, constraints):
            constrained_planets += 1.0/planet.probability

    percentage = (constrained_planets/total_planets) * 100.
    print constrained_planets, total_planets, percentage
    return percentage


def makeBarchart(weighted_occurences, N_bins, plot_range = [0,3000], tick_size=500, variable_name="T_planet", title = None):
    global FIGURE_COUNTER
    temp_list = weighted_occurences.keys()
    # Layout plot variables
    # font = {'size': 15}
    # plt.rc('font', **font)
    # font_size = '14'

    weighted_occurences_list = [0 for _ in range(N_bins)]
    planet_occurences_list = [0 for _ in range(N_bins)]
    starting_temp = plot_range[0]
    temp_range = plot_range[1] - plot_range[0]  #max(temp_list) - starting_temp
    temp_interval = temp_range / float(N_bins)

    HZ_count = 0
    pl_count = 0

    for temperature in temp_list:
        if plot_range[0] <= temperature <= plot_range[1]:
            weighted_occurences_list[int((temperature - starting_temp)/float(temp_interval))] += weighted_occurences[temperature]
            planet_occurences_list[int((temperature - starting_temp)/float(temp_interval))] += 1
            # pl_count += 1
            pl_count += weighted_occurences[temperature]
            if temperature >= 273 and temperature <= 373:
                HZ_count += weighted_occurences[temperature]
        else:
            continue

    x_labels = [starting_temp]
    x_positions = [0]
    temperature = starting_temp
    while temperature <= plot_range[1] - tick_size:  #max(temp_list):
        temperature += tick_size
        x_labels.append(temperature)
        x_positions.append((temperature - starting_temp )/ float(temp_interval))

    ### Plot percentages ###
    # for i in range(len(weighted_occurences_list)):
    #     weighted_occurences_list[i] = weighted_occurences_list[i]/float(pl_count)

    if title != None:
        plt.title(title)
    else:
        plt.title("Corrected Kepler planet distribution")

    # plt.xlabel("Temperature (K) ({0})".format(variable_name))
    plt.xlabel("Temperature (K)")
    plt.ylabel("Number of planets")
    plt.bar(left=range(N_bins), height=weighted_occurences_list, label="Nr of planets: " + str(int(round(pl_count)))) # "Planets in HZ : " + str(int(round(HZ_count))) + " of " + str(int(round(pl_count))))
    plt.legend(loc='upper right',prop={'size':13})
    plt.subplots_adjust(hspace=.4, wspace=.4)
    plt.xticks(x_positions, x_labels)
    return [starting_temp + i * temp_interval for i in range(N_bins + 1)], weighted_occurences_list



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
    x_list1 = []
    y_list1 = []
    x_list2 = []
    y_list2 = []
    x_list3 = []
    y_list3 = []
    x_list = []
    y_list = []
    total = 0
    HZ = 0
    HZ_incl_gas = 0
    for planet in planets:
        if eval("planet."+x_axis) == "" or eval("planet."+y_axis) == "" :
            continue
        if planet.M_planet == "":
            continue

        # try:
        #     x_list.append(float(eval("planet."+x_axis)))
        #     y_list.append(float(eval("planet."+y_axis)))
        # except:
        #     print (eval("planet."+x_axis))
        #     print "planet."+x_axis
        #     assert False
        total += 1

        if planet.phase == "liquid" and float(planet.R_planet) <= 2.5:
            x_list1.append(float(eval("planet."+x_axis)))
            y_list1.append(float(eval("planet."+y_axis)))
            HZ += 1
            HZ_incl_gas += 1
        # elif planet.phase != "":
        #     x_list2.append(float(eval("planet."+x_axis)))
        #     y_list2.append(float(eval("planet."+y_axis)))
        elif planet.phase == "liquid" and float(planet.R_planet) > 2.5:
            x_list2.append(float(eval("planet."+x_axis)))
            y_list2.append(float(eval("planet."+y_axis)))
            HZ_incl_gas += 1
        elif planet.phase == "gas":
            x_list2.append(float(eval("planet."+x_axis)))
            y_list2.append(float(eval("planet."+y_axis)))
        elif planet.phase == "solid":
            x_list3.append(float(eval("planet."+x_axis)))
            y_list3.append(float(eval("planet."+y_axis)))
        elif planet.phase == "supercritical fluid":
            x_list.append(float(eval("planet."+x_axis)))
            y_list.append(float(eval("planet."+y_axis)))


        # else:
        #     print planet.composition
        #     assert False

    plt.figure(FIGURE_COUNTER,  figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    FIGURE_COUNTER += 1
    font = {'size': 16}
    plt.rc('font', **font)
    plt.xlabel("Planet surface temperature (K)")
    plt.ylabel("$M_{planet}$ / $M_{Earth}$")
    x_range = np.linspace(0, 3000)
    y_range = np.linspace(0, 15)
    # plt.axvspan(273, 373, facecolor='g', alpha=0.5)
    plt.plot(x_range, [1 for i in x_range], c="k")
    plt.plot(x_range, [3.883 for i in x_range], c="k")
    plt.plot(x_range, [11.209 for i in x_range], c="k")
    # plt.plot([273 for i in y_range], y_range, "g--", lw=3)
    # plt.plot([373 for i in y_range], y_range, "g--", lw=3)
    plt.text(2980, 1.2,  'Earth', horizontalalignment='right')
    plt.text(2980, 3.883 + 0.2,  'Neptune', horizontalalignment='right')
    plt.text(2980, 11.209 + 0.2,  'Jupiter', horizontalalignment='right')
    # plt.xlabel(x_axis)
    # plt.ylabel(y_axis)
    if axis != None:
        plt.axis(axis)

    plt.xlim(0, 3000)
    plt.ylim(0, 20)
    # plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')

    # plt.scatter(x_list2, y_list2, color="#4682b4", label="Gas")
    # plt.scatter(x_list3, y_list3, color="#8b4513", label="Supercritical fluid")

    plt.scatter(x_list3, y_list3, color="black", label="Solid")

    plt.scatter(x_list2, y_list2, color="#4682b4", label="Gas")
    plt.scatter(x_list1, y_list1, color="green", label="Potential habitable planet")
    # plt.scatter(x_list, y_list, color="yellow", label="Supercritical fluid")
    plt.legend(bbox_to_anchor=(0.9, 0.84), bbox_transform=plt.gcf().transFigure)
    print 'Total planets = ', total
    print 'HZ planets = ', HZ, '       with gas planets = ', HZ_incl_gas
    print 'Percentage in HZ = ', (float(HZ)/float(total)) * 100
    print 'Percentage in HZ with gas planets = ', (float(HZ_incl_gas)/float(total)) * 100

def plotPropertyPercentage(planets, plot_range=(0,4,.3), variable = "R_star", HZ_extra = {"T_planet":(273, 373)}, constraints = {}):
    R_percentages = []
    x_list = []
    errors = []
    start, end, stepsize = plot_range
    number_points = []
    for i in np.arange(start + stepsize/2.0, end, stepsize):
        variable_dict = {variable:(float(i) - stepsize/2.0, float(i) + stepsize/2.0)}
        total_dict = dict(variable_dict.items() + constraints.items())
        total = getConstrainedPlanets(planets, total_dict)
        if len(total) < 5:
            print "not enough data points for {1} < {0} < {2}"\
                .format(variable, float(i) - stepsize/2.0, float(i) + stepsize/2.0)
            continue

        HZ_dict = dict(HZ_extra.items() + total_dict.items())
        HZ = getConstrainedPlanets(planets,  HZ_dict)
        for planet in HZ:
            assert planet in total

        values = []
        for planet in total:
            values.append(1.0 / planet.probability)

        perc = sum([1.0 / planet.probability for planet in HZ]) / sum([1.0 / planet.probability for planet in total]) * 100.
        R_percentages.append(perc)
        x_list.append(i)
        number_points.append((len(HZ), len(total)))
        errors.append(scipy.stats.sem(values, ddof=0))

    plt.scatter(x_list, R_percentages, s=100)
    plt.title("Habitable Zone percentage with {0} as variable".format(variable))
    font = {'size': 18}
    plt.rc('font', **font)
    plt.xlabel("Planet radius (R / R$_{Sun}$)")
    plt.ylabel("Percentage (%)")
    plt.errorbar(x_list, R_percentages, yerr=errors, linestyle="None") #, fmt="o")
    plt.ylim(0, max(R_percentages) + 20)
    # print calculatePercentage(data[0], {"T_planet":(273,373), "R_planet":(0.5,2)})

    for x,y,points in zip(x_list, R_percentages, number_points):
        plt.annotate(r"$\frac{%i}{%i}$"%(points[0], points[1]), (x + .03, y + 1))


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


def EarthMassToKg(mass):
    return float(mass) * 5.97219*10**24


def kgToEarthMass(mass):
    return mass /(5.97219*10**24)


def getDistance(mass, period):
    """
    Use Kepler's third law to calculate the distance between planet and star given M_star and P_orbit
    """
    return metersToAU((G * solarMasstoKg(mass) * (daysToSeconds(period)**2) / (4 * np.pi**2))**(1/3.))


def getPlanetTemperature(T_star_, R_star_, distance_, e=0, albedo=.3):
    """
    For a given star temperature and distance calculates the temperature (K) of the planet. Based on the equality of
    L_in and L_out.
    """

    R_star = float(solarRadiusToMeters(float(R_star_)))
    distance = float(AUToMeters(float(distance_)))
    T_star = float(T_star_)

    return ((1 - albedo)*(T_star**4)*(R_star**2)/(4 * (distance**2) * (1 - e)**(1/2.)))**(1/4.)

def makeStellarTypeSubplots(weighted_occurences_variable="T_planet", plott_range=[0,3000],tick=500):

    steller_type_list = ['A', 'K', 'M', 'G', 'F']
    plt.subplot(2,3,1)
    derived_data, derived_label_row  = ReadDerivedData()
    planets, combination, weighted_occurences = readData(derived_data, derived_label_row, weighted_occurences_variable)
    makeBarchart(weighted_occurences, BINS, plot_range=plott_range, tick_size=tick, variable_name="", title="All star types")

    for ster_type, i in zip(steller_type_list, range(2, len(steller_type_list) + 2)):
        plt.subplot(2, 3, i)
        derived_data, derived_label_row  = ReadDerivedData()
        planets, combination, weighted_occurences = readData(derived_data, derived_label_row, weighted_occurences_variable, constraints={"S_type": ster_type})
        makeBarchart(weighted_occurences, BINS, plot_range=plott_range, tick_size=tick, title= "" + ster_type + " Type star")


def safeDataMinitab():
    print "writing to file..."
    with open('minitlab_data.csv', 'wb') as csvfile:
            writer = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["T_planet", "M_planet", "M_star", "S_type", "distance", "probability", "weight"])
            for planet in planets:
                writer.writerow([planet.T_planet, planet.M_planet, planet.M_star, planet.S_type, planet.distance, planet.probability, gewicht])
    print "done"


def safeDataMinitab2(x_labels, heights):
    assert len(x_labels)  == len(heights)
    print "writing to file..."
    with open('minitlab_data.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["x", "height"])
        for i in range(len(heights)):
            height = heights[i]
            x_pos = x_labels[i]
            writer.writerow([x_pos, height])
    print "done"


def safeDataMathematica(planets):
    print "writing to file..."
    skips =0
    with open(r'C:\Users\joris\Dropbox\Workshop Sterrenkunde\HONGER.csv', 'wb') as csvfile:
            writer = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
            # writer.writerow(["ID_planet","T_planet", "P_planet"])
            for planet in planets:
                if planet.M_planet != "":
                    print planet.T_planet, planet.P_planet

                    writer.writerow([planet.ID, planet.T_planet, planet.P_planet])
                else:
                    skips += 1
    print "done but failed to write {0} planets".format(skips)


def gammaFit(x, theta=27967.73176, k=0.11456):
    # x /= 25.0
    return 1/(math.gamma(k) * theta**k) * x ** (k - 1) *\
           math.exp(-x/theta)


def gammaFit2(x, alpha, beta):
    return (beta**alpha) * (x ** (alpha - 1)) * np.exp(-x * beta) / (math.gamma(alpha))


def fit_gaussion(x, a, b, c):
    return a * np.exp((-1*(x - b)**2.)/(2*c**2.))


def fit_beta(x, a, b):
    return scipy.stats.betaprime.pdf(x, a, b)



derived_data, derived_label_row = readDerivedData()
phase_dict = readPhaseData()
planets, combination, weighted_occurences = readData(derived_data, derived_label_row, phase_dict)
fout = 0
fout_2 = 0
for planet in planets:
    if planet.phase == "liquid" and (373<planet.T_planet or planet.T_planet < 273):
        fout += 1
    elif planet.phase != "liquid" and (373>planet.T_planet> 273):
        fout_2 += 1

print abs(fout - fout_2)


# safeDataMathematica(planets)
# x, y = makeBarchart(weighted_occurences, BINS, plot_range=[0, 3000], tick_size=500, variable_name="T_planet")
# x_2 = [(x[i] + x[i+1])/2.0 for i in (range(len(x)-1))]
# plt.plot(x_2, y)
# plt.show()
# popt, pcov = scipy.optimize.curve_fit(fit_beta, np.array(x_2), np.array(y), p0=np.array([110, 200]))
# print popt, pcov
#
# plot_fit = []
# for i in x:
#     plot_fit.append(fit_beta(i, popt[0], popt[1]))


# plt.plot([val/BINS for val in x], plot_fit, c='r')
# makeBarchart(weighted_occurences, BINS, plot_range=[0, 3000], tick_size=500, variable_name="T_planet")

makeScatter(planets, "T_planet", "M_planet")
# plt.show()

# derived_data, derived_label_row  = ReadDerivedData()
# planets, combination, weighted_occurences = readData(derived_data, derived_label_row)
# x, y = makeBarchart(weighted_occurences, BINS, plot_range=[0,3000], tick_size=500, variable_name="T_planet")
# x_2 = [(x[i] + x[i+1])/2.0 for i in (range(len(x)-1))]
# fit_alpha,fit_loc,fit_beta=scipy.stats.gamma.fit(y, floc=0)
# # y = [gammaFit2(x, alpha=fit_alpha, beta=fit_beta) for x in x_2]
# y = [gammaFit(x, alpha=2, beta=6) for x in x_2]
# print y
# plt.figure(5)
# plt.plot(x_2,y)
# plt.show()


# y_3 = []
# for x in np.arange(0,20,.001):
#     y_3.append(gammaFit(x))
# plt.plot(np.arange(0,20,.001), y_3)



# makeStellarTypeSubplots()


# plt.figure(1)
# derived_data, derived_label_row  = ReadDerivedData()
# planets, combination, weighted_occurences = readData(derived_data, derived_label_row)
# makeBarchart(weighted_occurences, BINS)
#
# plt.figure(0)
# derived_data, derived_label_row  = ReadDerivedData()
# planets, combination, weighted_occurences = readData(derived_data, derived_label_row, weighted_occurences_variable="R_planet")
# makeBarchart(weighted_occurences, BINS, plot_range=[0,15], tick_size=3, variable_name="R_planet")
#
# plt.figure(2)

# derived_data, derived_label_row  = ReadDerivedData()
# planets, combination, weighted_occurences = readData(derived_data, derived_label_row, weighted_occurences_variable="density", constraints={"M_planet": (0, 10**5)})
# makeBarchart(weighted_occurences, BINS, plot_range=[1000,10000], tick_size=1500, variable_name="density")
#


# makeBarchart(data[2], BINS, [0,500], 100)
# makeScatter(data[0])
# makeBarchart(data[2], BINS)
# makeScatter(data[0])
# makeScatter(data[0], axis=[0,3000, 0, 10])
# makeHistogram(data[0], data[1])
plt.show()

# derived_data, derived_label_row  = ReadDerivedData()
# planets, combination, weighted_occurences = readData(derived_data, derived_label_row, weighted_occurences_variable="R_planet")
# plotPropertyPercentage(planets, variable="R_planet", plot_range=(0,5,.8))
# # def plotPropertyPercentage(planets, plot_range=(0,4,.3), variable = "R_star", HZ_extra = {"T_planet":(273, 373)}, constraints = {}):
# plt.show()
