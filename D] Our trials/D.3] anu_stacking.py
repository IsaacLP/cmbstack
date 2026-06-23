

# These are our toolboxes
import healpy as hp   # This handles the round sky maps
import numpy as np    # This handles numbers and math
import matplotlib.pyplot as plt  # This draws pictures

# ----------------------------------------------------------

# CMB Peak Stacking Pipeline

## Overview
#This repository contains our Python implementation for stacking patches of the Cosmic Microwave Background (CMB) temperature sky. Our goal is to take a theoretical power spectrum, generate a synthetic sky map, locate the hottest local maxima (peaks), and average the surrounding patches. This averaging process, known as "stacking," enhances the symmetric signal of the typical peak profile while suppressing random noise and cosmic variance.

#Below, we document the mathematical reasoning behind each of our six core functions.

## Summary of Workflow
#1. **Input:** Theoretical $D_\ell$ spectrum.
#2. **Convert:** Change $D_\ell$ to $C_\ell$.
#3. **Simulate:** Generate a random Gaussian CMB sky map.
#4. **Standardize:** Remove the mean and normalize by $\sigma$.
#5. **Identify:** Find and rank the local maxima.
#6. **Extract:** Cut out $5^\circ$ discs around the top 100 peaks.
#7. **Average:** Stack the patches to reveal the universal peak profile.
#—----------------------------------------------------------------------------------------------


# -------------------------------------------------------------

# FUNCTION 1: Read the data file from your supervisor

### 1. Reading the Power Spectrum (Converting $D_\ell$ to $C_\ell$)
# **Function:** `read_power_spectrum_from_file`

#The data provided by our supervisors contains the power spectrum in the form of $D_\ell^{TT}$, defined as:
#\[
#D_\ell \equiv \frac{\ell(\ell+1)}{2\pi} C_\ell
#\]
#However, the Healpy library expects the angular power spectrum $C_\ell$ as input for generating maps. Therefore, we rearrange the formula to convert our data:
#\[
#C_\ell = \frac{2\pi}{\ell(\ell+1)} D_\ell
#\]
#We set $C_0 = 0$ to avoid the division-by-zero singularity at $\ell=0$, as the monopole contributes no fluctuating signal.




# ----------------------------------------------------------
def read_power_spectrum_from_file(filename):
    """
    Opens the file, skips the header line, 
    and grabs the Temperature (TT) column.
    """
    # Load the file, skip the first line (which has the # and headers)
    raw_data = np.loadtxt(filename, skiprows=1)
    
    # Column 0 = the scale of the bumps
    scale_of_bumps = raw_data[:, 0]
    
    # Column 1 = the Temperature power values (this is the only one we care about)
    temperature_power_values = raw_data[:, 1]
    
    # Convert to the format Healpy wants (don't worry about the math, just do it)
    # Formula: needed_format = temperature_power_values * 2 * pi / (scale * (scale+1))
    temperature_power_spectrum = temperature_power_values * 2 * np.pi / (scale_of_bumps * (scale_of_bumps + 1))
    
    # The first one causes a division-by-zero, so we just set it to 0
    temperature_power_spectrum[0] = 0
    
    return temperature_power_spectrum


# ----------------------------------------------------------
# FUNCTION 2: Create a fake CMB sky map

### 2. Generating the Synthetic CMB Map
#**Function:** `create_fake_sky_map`

#We generate a random, Gaussian-distributed realization of the CMB sky. Healpy achieves this by drawing spherical harmonic coefficients $a_{\ell m}$ from a Gaussian distribution with variance $C_\ell$, formally defined as:
#\[
#\langle |a_{\ell m}|^2 \rangle = C_\ell
#\]
#The temperature at any point on the sphere $\hat{n}$ is then computed via the spherical harmonic transform:
#\[
#T(\hat{n}) = \sum_{\ell=0}^{\ell_{\text{max}}} \sum_{m=-\ell}^{\ell} a_{\ell m} \, Y_{\ell m}(\hat{n})
#\]
#where $Y_{\ell m}$ are the spherical harmonic functions. The `synfast` routine performs this exact computation to produce our mock sky map.


# ----------------------------------------------------------
def create_fake_sky_map(temperature_power_spectrum, sky_resolution=128):
    """
    Makes a picture of the whole sky using the power spectrum.
    sky_resolution=128 means a decent quality picture.
    """
    fake_sky = hp.synfast(temperature_power_spectrum, sky_resolution)
    return fake_sky


# ----------------------------------------------------------
# FUNCTION 3: Clean up the map (remove average & scale it)

### 3. Preprocessing the Map (Removing Monopole & Normalizing)
#**Function:** `remove_average_and_scale_map`

#Before searching for peaks, we must standardize our map to ensure the detection threshold is uniform. First, we remove the monopole (the average temperature) from the map:
#\[
#T_{\text{fluctuation}}(\hat{n}) = T(\hat{n}) - \langle T \rangle
#\]
#Next, we divide by the standard deviation $\sigma$ to make the map dimensionless and unit-variant:
#\[
#T_{\text{norm}}(\hat{n}) = \frac{T_{\text{fluctuation}}(\hat{n})}{\sigma}
#\]
#where $\sigma = \sqrt{\langle (T - \langle T \rangle)^2 \rangle}$. After this step, our map has a mean of zero and a standard deviation of one, meaning peaks are measured in units of "sigma" above the average.



# ----------------------------------------------------------
def remove_average_and_scale_map(sky_map):
    """
    Step 1: Subtract the average temperature (so the map centers around 0).
    Step 2: Divide by the standard deviation (so the map has nice, small numbers).
    """
    sky_minus_average = sky_map - np.mean(sky_map)
    sky_normalized = sky_minus_average / np.std(sky_minus_average)
    return sky_normalized


# ----------------------------------------------------------
# FUNCTION 4: Find the hottest spots (peaks)

### 4. Detecting the Hottest Local Maxima
#**Function:** `find_hottest_spots`

#We utilize Healpy's `hotspots` algorithm to identify local maxima. A pixel is considered a local maximum if its temperature value is strictly greater than all its immediate neighboring pixels in the HEALPix grid.

#Once all maxima are found, we rank them by their normalized temperature $T_{\text{norm}}$. We select the top $N=100$ hottest spots:
#\[
#\text{Peak Set} = \{ \hat{n}_p \in \text{Maxima} \mid T_{\text{norm}}(\hat{n}_p) \text{ is in the top } N \}
#\]
#We focus on the hottest peaks because they represent the most significant structures in the Gaussian field and yield the highest signal-to-noise ratio in the final stack.


# ----------------------------------------------------------
def find_hottest_spots(cleaned_sky_map, number_of_peaks_to_find=100):
    """
    Finds all the local hot spots, picks the brightest ones, 
    and returns their pixel locations.
    """
    # This finds every single hot spot on the map
    everything, cold_pixels, all_hot_pixels = hp.hotspots(cleaned_sky_map)
    
    # Check the temperature at each hot spot
    temperatures_at_hot_spots = cleaned_sky_map[all_hot_pixels]
    
    # Sort them so the hottest comes first
    order_from_hottest_to_coldest = np.argsort(temperatures_at_hot_spots)[::-1]
    
    # Keep only the top ones (e.g., top 100)
    hottest_peak_pixels = all_hot_pixels[order_from_hottest_to_coldest][:number_of_peaks_to_find]
    
    return hottest_peak_pixels


# ----------------------------------------------------------
# FUNCTION 5: Cut out little circles around those spots

### 5. Extracting Patches (Cutting Circles Around Peaks)
#**Function:** `cut_circles_around_spots`

#For each selected peak located at direction vector $\hat{n}_p$, we define a circular aperture (patch) with a fixed angular radius $\Theta = 5^\circ$. 

#To extract the patch, we find every pixel direction $\hat{n}$ on the sphere that falls within this radius. The angular distance $\Delta \theta$ between the peak center and a given pixel is computed via the dot product:
#\[
#\Delta \theta = \arccos\left( \hat{n}_p \cdot \hat{n} \right)
#\]
#We include all pixels satisfying $\Delta \theta \leq \Theta$. This gives us a set of temperature values $T_{\text{norm}}(\hat{n})$ surrounding each peak. Conceptually, we are recentering these patches so that the peak lies exactly at the origin of our local coordinate system.


# ----------------------------------------------------------
def cut_circles_around_spots(cleaned_sky_map, hottest_peak_pixels, sky_resolution, circle_size_in_degrees=5):
    """
    For each hot spot, draws a circle (like a cookie cutter) 
    and saves all the temperatures inside that circle.
    """
    # Convert degrees to radians (because Healpy speaks radians)
    circle_size_in_radians = np.radians(circle_size_in_degrees)
    
    all_circles = []  # This will hold all the little cut-out circles
    
    # Loop through each hot spot
    for pixel_location in hottest_peak_pixels:
        # Get the latitude and longitude of this spot
        latitude, longitude = hp.pix2ang(sky_resolution, pixel_location)
        
        # Convert that to a 3D direction arrow (x,y,z)
        center_direction_vector = hp.ang2vec(latitude, longitude)
        
        # Find all the pixel-numbers that fall inside this circle
        pixels_in_circle = hp.query_disc(sky_resolution, center_direction_vector, circle_size_in_radians, inclusive=True)
        
        # Get the actual temperature values at those pixels
        temperatures_in_circle = cleaned_sky_map[pixels_in_circle]
        
        # Save this circle to our list
        all_circles.append(temperatures_in_circle)
    
    return all_circles


# ----------------------------------------------------------
# FUNCTION 6: Average all those circles together (STACKING!)


### 6. Stacking (Averaging the Patches)
#**Function:** `average_all_circles`

#This is the core of our project. To construct the mean stacked peak map, we average the extracted patches pixel-by-pixel. For $N$ peaks, the stacked radial profile is defined as:
#\[
#S(\Delta \theta) = \frac{1}{N} \sum_{i=1}^{N} T_{\text{norm}}^{(i)}(\Delta \theta)
#\]
#where $i$ indexes the selected peaks. 

#By averaging over many independent patches, random fluctuations and noise (which are uncorrelated) average down towards zero, while the coherent, symmetric profile of a typical CMB hot spot remains. This final average gives us the characteristic "peak shape" of the CMB temperature anisotropies.
# ----------------------------------------------------------
def average_all_circles(all_circles):
    """
    Takes all the circles, makes them the exact same size, 
    and averages them pixel-by-pixel to create ONE "stacked" circle.
    """
    # Find the smallest circle (just to make sure they all fit together)
    smallest_circle_size = min(len(circle) for circle in all_circles)
    
    # Trim every circle to that exact size
    circles_same_size = [circle[:smallest_circle_size] for circle in all_circles]
    
    # Turn the list into a table (rows = circles, columns = pixels)
    table_of_circles = np.array(circles_same_size)
    
    # Average down the rows (axis=0 means "average each column across all rows")
    average_circle = np.mean(table_of_circles, axis=0)
    
    return average_circle


# ==========================================================
# THIS IS THE MAIN PROGRAM - IT RUNS EVERYTHING IN ORDER
# ==========================================================

if __name__ == "__main__":
    
    print("Step 1: Reading your data file...")
    # !!! CHANGE THIS TO YOUR ACTUAL FILENAME !!!
    my_data_file = "copy_base_plikHM_TTTEEE_lowl_lowE_lensing.minimum.txt"  
    power_spectrum = read_power_spectrum_from_file(my_data_file)

    print("Step 2: Generating a fake sky map...")
    my_sky_map = create_fake_sky_map(power_spectrum)

    print("Step 3: Cleaning up the map (removing average and scaling)...")
    cleaned_sky = remove_average_and_scale_map(my_sky_map)

    print("Step 4: Finding the 100 hottest spots on the sky...")
    hot_spots = find_hottest_spots(cleaned_sky, number_of_peaks_to_find=100)
    print(f"Found {len(hot_spots)} hot spots!")

    print("Step 5: Cutting out circles around those hot spots...")
    my_circles = cut_circles_around_spots(cleaned_sky, hot_spots, sky_resolution=128)

    print("Step 6: Averaging all the circles together (this is the STACK!)...")
    stacked_result = average_all_circles(my_circles)

    print("All done! Drawing a picture of the stacked result...")
    
    # Draw a simple chart showing the temperatures inside the stacked circle
    plt.figure(figsize=(6, 5))
    plt.hist(stacked_result, bins=30, edgecolor='black', color='skyblue')
    plt.title("The Stacked Peak! (Average temperature of all hot spots)")
    plt.xlabel("Temperature value (normalized)")
    plt.ylabel("Number of pixels")
    plt.grid(True, alpha=0.3)
    plt.show()
    
    print("Check the pop-up window for your plot. Success!")
