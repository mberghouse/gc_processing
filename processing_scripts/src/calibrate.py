import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt


# Function to extract peak table info from Shimadzu
# GC-2014 ascii text output
def get_peak_info(filename, N2O_check=1, CH4_check=1, CO2_check=1):
    N2O_area = 0
    CO2_area = 0
    CH4_area = 0

    # Create check to make sure (1 = peak was found in file)
    [N2O_found, CO2_found, CH4_found] = [0, 0, 0]

    f = open(filename)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) > 9:
            # print(line)
            if fields[10] == 'Nitrous oxide':
                N2O_area = int(fields[4])
                N2O_found = 1
            elif fields[10] == 'CO2':
                CO2_area = int(fields[4])
                CO2_found = 1
            elif fields[10] == 'Methane':
                CH4_area = int(fields[4])
                N2O_found = 1

    f.close()

    # Print out filename and unfound peaks (do not return error)
    if N2O_found != 1 and N2O_check == 1:
        print(filename+": N2O peak not found!")
    if CH4_found != 1 and CH4_check == 1:
        print(filename+": CH4 peak not found!")
    if CO2_found != 1 and CO2_check == 1:
        print(filename+": CO2 peak not found!")

    return N2O_area, CO2_area, CH4_area


# Given a dataframe with GC filenames as a single col (in the appropriate format),
# return a copy of that dataframe with the location, date, depth and vial number
# in separate columns
def sample_info_from_filename(df, column_name):
    # Create dataframe of all samples
    samples = df[df[column_name].str.match('Gas')].copy()

    # Get the vial position in the autoloader
    samples['vial_pos'] = samples[column_name].str.extract(r'(\d+)\.t') # Extract digits before '.t', not including '.t'

    # Get sampling method
    samples['method'] = samples[column_name].str.extract(r'-([GH])_')

    # Get the sampling date
    samples['date'] = samples[column_name].str.extract(r'(201.....)')

    # Get location
    samples['location'] = samples[column_name].str.extract(r'Gas_([A-Za-z]+[a-z0-9])-[0-9]{1,3}cm') # Extract up to 'cm'

    # Get depth
    samples['depth'] = samples[column_name].str.extract(r'-([0-9]{1,3})cm')

    return samples


# Define a function to calculate m, b and r-squared for
# a simple linear regression
def linreg(x, y):
    regres = {}
    coeffs = np.polyfit(x, y, 1)

    regres['poly'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)  # or [p(z) for z in x]
    ybar = np.sum(y) / len(y)  # or sum(y)/len(y)
    ssreg = np.sum((yhat - ybar) ** 2)  # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar) ** 2)  # or sum([ (yi - ybar)**2 for yi in y])
    regres['r-squared'] = ssreg / sstot

    return regres


# Define a quadratic fitting function
def quadreg(x, y):
    regres = {}
    coeffs = np.polyfit(x, y, 2)
    regres['poly'] = coeffs.tolist()

    # Dummy value for r-squared
    regres['r-squared'] = float('nan')

    return regres


def quadratic(xinput, coefficients):
    return coefficients[0] * xinput**2 + coefficients[1] * xinput + coefficients[2]


# Plot the standards showing the regression trendline and
# the r-squared value on the plot
def plot_standards(stdx, stdy, label, power=1):
    # Apply a simple least squares regression
    if power == 1:
        regres = linreg(stdx, stdy)
        p = np.poly1d(regres['poly'])
        xhat = stdx
        yhat = p(xhat)
    elif power == 2:
        regres = quadreg(stdx, stdy)
        xhat = np.linspace(min(stdx), max(stdx), num=50)
        yhat = quadratic(xhat, regres['poly'])
    else:
        print("ERROR --> Power must equal 1 or 2")
        return -1

    plt.figure(figsize=(16, 8))
    plt.scatter(stdx, stdy, edgecolors='r', facecolors='none', marker='s', s=120, label="Std")
    plt.plot(xhat, yhat, "k-")

    # Add correlation coefficient to the plot
    x_loc = pd.Series.min(stdx)
    y_loc = pd.Series.max(stdy) * .9
    rform = "{:1.4f}".format(regres['r-squared'])
    plt.text(x_loc, y_loc, r'$R^2$ = ' + rform)

    # Add axis labels and title
    plt.xlabel(label + ' Concentration (ppm)')
    plt.ylabel(r'Peak area')
    plt.title(label)
    plt.legend(loc='best')


# Plot the standards along with the samples
def plot_samples(x, stdx, stdy, label, power=1):

    if power == 1:
        # Apply a simple linear regression to the standards
        regres = linreg(stdx, stdy)
        p = np.poly1d(regres['poly'])
        m, b = regres['poly']
        xhat = stdx
        yhat = p(xhat)
        predx = (x - b) / m
    elif power == 2:
        regres = quadreg(stdx, stdy)
        qcoef = regres['poly'].copy()
        xhat = np.linspace(min(stdx), max(stdx), num=50)
        yhat = quadratic(xhat, qcoef)

        # Make a prediction by finding the roots of the adjusted polynomial
        predx = x.tolist()
        i = 0
        for speak in x:
            scoef = qcoef.copy()
            scoef[2] = scoef[2] - speak
            roots = np.roots(scoef)
            predx[i] = roots[roots < 80]
            i += 1
    else:
        print("ERROR --> Power must equal 1 or 2")
        return -1

    # Plot the standards
    plt.figure(figsize=(16, 8))
    plt.scatter(stdx, stdy, edgecolors='r', facecolors='none', marker='s', s=120, label="Std")
    plt.plot(xhat, yhat, "k-")

    # Plot the samples
    plt.scatter(predx, x, edgecolors='b', facecolors='b', marker='o', s=120, label="Samples")

    # Add axis labels and title
    plt.xlabel(label + ' Concentration (ppm)')
    plt.ylabel(r'Peak area')
    plt.title(label)
    plt.legend(loc='best')


# Since the data look good, export the sample concentrations to a text file
# File format:
# [NUMBER] [LOCATION] [DEPTH] [DATE] [N2O] [CO2]
def save_sampleconc(N2Ostdx, N2Ostdy, CO2stdx, CO2stdy, CH4stdx, CH4stdy,
                    samples, save_path, write_to_datafile=False):
    # Apply a simple least squares regression to the standards
    Nregres = linreg(N2Ostdx, N2Ostdy)
    Nm, Nb = Nregres['poly']
    N2Osam = samples.loc[:, "N2O_peak"].astype(float)

    Cregres = linreg(CO2stdx, CO2stdy)
    Cm, Cb = Cregres['poly']
    CO2sam = samples.loc[:, "CO2_peak"].astype(float)

    if len(CH4stdx) > 0:
        Hregres = linreg(CH4stdx, CH4stdy)
        Hm, Hb = Hregres['poly']
        CH4sam = samples.loc[:, "CH4_peak"].astype(float)

    # Create a new dataframe in which to accumulate the output
    # First, check if the sample name is formatted correctly
    pattern = re.compile('[0-9]{3}_NT[T,C][0-9]_[0-9]{3}cm_[0-9]{8}')
    for name in samples['Sample_name']:
        if not re.match(pattern, name):
            print("Warning - This sample name is not formatted correctly:")
            print(name)
            return

    output = pd.DataFrame(samples['Sample_name'].str.slice(0, 3))  # Sample number
    output.columns = ['SampleID']
    output['location'] = pd.Series(samples['Sample_name'].str.slice(4, 8))  # Sample location
    output['depth'] = pd.Series(samples['Sample_name'].str.slice(9, 12).astype(float))  # Convert depth to float
    output['date'] = pd.Series(
        pd.to_datetime(samples['Sample_name'].str.slice(15, 23)))  # Convert date to datetime
    output['dup'] = pd.Series(pd.Series(samples['Sample_name'].str.contains('dup')))  # Bool of duplicate or not
    output['method'] = 'G'
    output['vial'] = 0

    # Set detection limits to lowest standard value
    if len(CH4stdx) > 0:
        CH4dl = min(CH4stdx)
    if len(N2Ostdx) > 0:
        N2Odl = min(N2Ostdx)
    if len(CO2stdx) > 0:
        CO2dl = min(CO2stdx)

    # Calculate N2O conc if we have standards, otherwise nan
    if len(N2Ostdx) > 0:
        output['N2O_conc'] = pd.Series((N2Osam - Nb) / Nm)
        # Set any values below the detection limit to 1/2 the detection limit
        output['N2O_conc'] = output['N2O_conc'].where(output['N2O_conc'] > N2Odl, N2Odl*0.5)
        output['N2O_bdl'] = output['N2O_conc'] <= N2Odl
    else:
        output['N2O_conc'] = np.nan
        output['N2O_bdl'] = np.nan

    # Calculate CO2 conc if we have standards, otherwise nan
    if len(CO2stdx) > 0:
        output['CO2_conc'] = pd.Series((CO2sam - Cb) / Cm)
        # Set any values below the detection limit to 1/2 the detection limit
        output['CO2_conc'] = output['CO2_conc'].where(output['CO2_conc'] > CO2dl, CO2dl*0.5)
        output['CO2_bdl'] = output['CO2_conc'] <= CO2dl
    else:
        output['CO2_conc'] = np.nan
        output['CO2_bdl'] = np.nan

    # Calculate CH4 conc if we have standards, otherwise nan
    if len(CH4stdx) > 0:
        output['CH4_conc'] = pd.Series((CH4sam - Hb) / Hm)
        # Set any values below the detection limit to 1/2 the detection limit
        output['CH4_conc'] = output['CH4_conc'].where(output['CH4_conc'] > CH4dl,  CH4dl*0.5)
        output['CH4_bdl'] = output['CH4_conc'] <= CH4dl
    else:
        output['CH4_conc'] = np.nan
        output['CH4_bdl'] = np.nan

    # Rearrange column order
    cols = ['vial', 'location', 'depth', 'date', 'method', 'dup', 'N2O_conc', 'N2O_bdl',
                                                                  'CO2_conc', 'CO2_bdl',
                                                                  'CH4_conc', 'CH4_bdl']
    output = output[cols]

    # Append the output to a text file if requested
    if write_to_datafile:
        output.to_csv(save_path, mode='w', sep="\t", index=False, header=False, na_rep='NaN')

    return output


# Since the data look good, export the sample concentrations to a text file
# File format:
# [NUMBER] [LOCATION] [DEPTH] [DATE] [N2O] [CO2]
def save_sampleconc_new(N2Ostdx, N2Ostdy, CO2stdx, CO2stdy, CH4stdx, CH4stdy,
                    samples, save_path, write_to_datafile=False):
    # Apply a simple least squares regression to the standards
    if len(N2Ostdx) > 0:
        Nregres = linreg(N2Ostdx, N2Ostdy)
        Nm, Nb = Nregres['poly']
        N2Osam = samples.loc[:, "N2O_peak"].astype(float)

    if len(CO2stdx) > 0:
        Cregres = linreg(CO2stdx, CO2stdy)
        Cm, Cb = Cregres['poly']
        CO2sam = samples.loc[:, "CO2_peak"].astype(float)

    if len(CH4stdx) > 0:
        Hregres = linreg(CH4stdx, CH4stdy)
        Hm, Hb = Hregres['poly']
        CH4sam = samples.loc[:, "CH4_peak"].astype(float)

    # Create a new dataframe in which to accumulate the output
    # First, check if the sample name is formatted correctly
    pattern = re.compile('Gas_[A-Za-z0-9]+-[0-9]+cm_[0-9]{8}-[H,G]_[0-9]+_[0-9]+.txt')
    for name in samples['sample']:
        if not re.match(pattern, name):
            print("Warning - This sample name is not formatted correctly:")
            print(name)
            return

    # Split name by underscore
    output = samples['sample'].str.split('_', expand=True)
    output.columns = ['_', 'sampler', 'dateandmethod', 'run_date', 'vial']

    # Split location/depth by dash
    tmp = output['sampler'].str.split('-', expand=True)
    output['location'] = tmp[0]
    output['depth'] = tmp[1]
    output['depth'] = output['depth'].str.slice(0, -2)  # remove trailing 'cm'

    # Split date/method by dash
    tmp = output['dateandmethod'].str.split('-', expand=True)
    output['date'] = pd.to_datetime(tmp[0])
    output['method'] = tmp[1]
    output['vial'] = output['vial'].str.slice(0, -4)  # remove trailing '.txt'
    output['dup'] = samples['sample'].str.contains('dup')  # Bool of duplicate or not

    # Drop extraneous columns
    output.drop(columns=['_', 'sampler', 'dateandmethod'], inplace=True)

    # Set detection limits to lowest standard value
    if len(CH4stdx) > 0:
        CH4dl = min(CH4stdx)
    if len(N2Ostdx) > 0:
        N2Odl = min(N2Ostdx)
    if len(CO2stdx) > 0:
        CO2dl = min(CO2stdx)

    # Now calculate concentrations
    # Calculate N2O conc if we have standards, otherwise nan
    if len(N2Ostdx) > 0:
        output['N2O_conc'] = pd.Series((N2Osam - Nb) / Nm)
        # Set any values below the detection limit to 1/2 the detection limit
        output['N2O_conc'] = output['N2O_conc'].where(output['N2O_conc'] > N2Odl, N2Odl*0.5)
        output['N2O_bdl'] = output['N2O_conc'] <= N2Odl
    else:
        output['N2O_conc'] = np.nan
        output['N2O_bdl'] = np.nan

    # Calculate CO2 conc if we have standards, otherwise nan
    if len(CO2stdx) > 0:
        output['CO2_conc'] = pd.Series((CO2sam - Cb) / Cm)
        # Set any values below the detection limit to 1/2 the detection limit
        output['CO2_conc'] = output['CO2_conc'].where(output['CO2_conc'] > CO2dl, CO2dl*0.5)
        output['CO2_bdl'] = output['CO2_conc'] <= CO2dl
    else:
        output['CO2_conc'] = np.nan
        output['CO2_bdl'] = np.nan

    # Calculate CH4 conc if we have standards, otherwise nan
    if len(CH4stdx) > 0:
        output['CH4_conc'] = pd.Series((CH4sam - Hb) / Hm)
        # Set any values below the detection limit to 1/2 the detection limit
        output['CH4_conc'] = output['CH4_conc'].where(output['CH4_conc'] > CH4dl,  CH4dl*0.5)
        output['CH4_bdl'] = output['CH4_conc'] <= CH4dl
    else:
        output['CH4_conc'] = np.nan
        output['CH4_bdl'] = np.nan

    # Rearrange column order
    cols = ['vial', 'location', 'depth', 'date', 'method', 'dup', 'N2O_conc', 'N2O_bdl',
                                                                  'CO2_conc', 'CO2_bdl',
                                                                  'CH4_conc', 'CH4_bdl']
    output = output[cols]

    # Append the output to a text file if requested
    if write_to_datafile:
        output.to_csv(save_path, mode='a', sep="\t", index=False, header=False, na_rep='NaN')

    return output
