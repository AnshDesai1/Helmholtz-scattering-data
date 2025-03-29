import numpy as np

# Image functions as matrices
def is_inside_circle(y1, y2, radius, y1c, y2c):
    # Shift coordinates by the circle center
    y1_shifted = y1 - y1c
    y2_shifted = y2 - y2c
    # Check circle equation
    return (y1_shifted**2 + y2_shifted**2) <= radius**2

# Matrix Form of Circle
def circle_matrix(xlim,Ngrid,numcircle,circle_data,nval):
    y1 = np.linspace(-xlim, xlim, Ngrid)
    y2 = np.linspace(-xlim, xlim, Ngrid)
    Y1, Y2 = np.meshgrid(y1,y2)
    image = np.zeros((Ngrid,Ngrid), dtype=float)
    for i in range(1, numcircle+1):
        R = circle_data[f'R{i}']
        y1c, y2c = circle_data[f'xcen{i}']
        value = nval[i - 1] - 1
        inside = is_inside_circle(Y1,Y2,R,y1c,y2c)
        image = np.where(inside,value,image)
    return image
                                  
def is_inside_ellipse(y1, y2, a, b, y1c, y2c, angle):
    y1_shifted = y1 - y1c # Shift horizontal domain
    y2_shifted = y2 - y2c # Shift vertical domain
    cos_angle = np.cos(-angle)
    sin_angle = np.sin(-angle)
    y1_rotated = y1_shifted * cos_angle + y2_shifted * sin_angle # Rotate coordinates
    y2_rotated = -y1_shifted * sin_angle + y2_shifted * cos_angle
    return (y1_rotated**2 / a**2 + y2_rotated**2 / b**2) <= 1
    
def ellipses_matrix(xlim, Ngrid, numellip, ellip_data, nval):
    y1 = np.linspace(-xlim, xlim, Ngrid)  # Discretize horizontal domain
    y2 = np.linspace(-xlim, xlim, Ngrid)  # Discretize vertical domain
    Y1, Y2 = np.meshgrid(y1, y2)  # Produce matrices of coordinates on square domain
    image = np.zeros((Ngrid, Ngrid), dtype=float)
    for i in range(1, numellip + 1):
        a = ellip_data[f'R{i}a']
        b = ellip_data[f'R{i}b']
        y1c, y2c = ellip_data[f'xcen{i}']
        angle = ellip_data[f'ang{i}']
        value = nval[i - 1] - 1
        inside = is_inside_ellipse(Y1, Y2, a, b, y1c, y2c, angle)
        image = np.where(inside, value, image)
    return image
    
def ellipses_circ_matrix(xlim, Ngrid, numellip, ellip_data, nval, outer_nval, outer_rad):
    y1 = np.linspace(-xlim, xlim, Ngrid)  # Discretize horizontal domain
    y2 = np.linspace(-xlim, xlim, Ngrid)  # Discretize vertical domain
    Y1, Y2 = np.meshgrid(y1, y2)  # Produce matrices of coordinates on square domain
    # Initialize image with the value of the outer circle
    image = np.ones((Ngrid, Ngrid), dtype=float) * (outer_nval - 1)
    # Set values inside the outer circle to outer_nval - 1, and outside to 0
    inside_circle = (Y1**2 + Y2**2) <= outer_rad**2
    image = np.where(inside_circle, outer_nval - 1, 0)
    # Iterate over ellipses
    for i in range(1, numellip + 1):
        print(f"{i} of {numellip}")
        a = ellip_data[f'R{i}a']
        b = ellip_data[f'R{i}b']
        y1c, y2c = ellip_data[f'xcen{i}']
        angle = ellip_data[f'ang{i}']
        value = nval[i - 1] - 1
        # Check if points are inside the ellipse
        inside = is_inside_ellipse(Y1, Y2, a, b, y1c, y2c, angle)
        # Update image only where the current ellipse is inside and hasn't been assigned yet
        image = np.where(np.logical_and(inside, image == (outer_nval - 1)), value, image)
    return image

def blob_matrix(Ngrid, xlim, Rmax, Rmin, npoints, angs, mags, nval):    
    Rads = (Rmax - Rmin) * mags + Rmin
    xav = np.mean(Rads * np.cos(angs))
    yav = np.mean(Rads * np.sin(angs))
    x_points = []
    y_points = []
    for j in range(npoints):
        x = Rads[j] * np.cos(angs[j]) - xav
        y = Rads[j] * np.sin(angs[j]) - yav
        x_points.append(x)
        y_points.append(y)
    blob = np.vstack((x_points, y_points)).T
    y1 = np.linspace(-xlim, xlim, Ngrid)
    y2 = np.linspace(-xlim, xlim, Ngrid)
    Y1, Y2 = np.meshgrid(y1, y2)
    grid_points = np.vstack((Y1.flatten(), Y2.flatten())).T
    blob_path = Path(blob)
    inside_blob = blob_path.contains_points(grid_points)
    image = np.zeros((Ngrid, Ngrid))
    image.flat[inside_blob] = nval - 1
    return image

def U_matrix(Box, xcen, t, ang, Ngrid, nval, xlim):
    # Unpack parameters
    Z, W, X, Y = Box
    a, b = xcen
    # Define grid extent based on the bounding box
    x = np.linspace(-1, 1, Ngrid)
    y = np.linspace(-1, 1, Ngrid)
    grid_x, grid_y = np.meshgrid(x, y)
    # Translate grid to center (a, b)
    grid_x -= a
    grid_y -= b
    # Apply inverse rotation to grid
    cos_ang = np.cos(ang+np.pi/2)
    sin_ang = np.sin(ang+np.pi/2)
    x_rot = cos_ang * grid_x - sin_ang * grid_y
    y_rot = sin_ang * grid_x + cos_ang * grid_y
    # Initialize the matrix
    matrix = np.zeros((Ngrid, Ngrid), dtype=float)
    # Define the U-shape in the rotated frame
    # Left vertical bar of U
    left_bar = (X <= x_rot) & (x_rot <= X + t) & (Z <= y_rot) & (y_rot <= W)
    # Right vertical bar of U
    right_bar = (Y - t <= x_rot) & (x_rot <= Y) & (Z <= y_rot) & (y_rot <= W)
    # Bottom horizontal bar of U
    bottom_bar = (X <= x_rot) & (x_rot <= Y) & (Z <= y_rot) & (y_rot <= Z + t)
    # Debug: Print number of points in each part
    print(f"Left bar points: {np.sum(left_bar)}")
    print(f"Right bar points: {np.sum(right_bar)}")
    print(f"Bottom bar points: {np.sum(bottom_bar)}")
    # Combine all parts of the U-shape
    u_shape = left_bar | right_bar | bottom_bar
    # Remove the hollow center
    hollow_center = (X + t <= x_rot) & (x_rot <= Y - t) & (Z + t <= y_rot) & (y_rot <= W - t)
    u_shape = u_shape & ~hollow_center
    # Assign the value nval - 1 to the U-shape
    matrix[u_shape] = nval - 1
    return matrix

def shepp_matrix(Ngrid, data, nval):
    matrix = np.zeros((Ngrid, Ngrid))  # Initialize background to 0
    # Create a grid from -1 to 1
    x = np.linspace(-1, 1, Ngrid)
    y = np.linspace(-1, 1, Ngrid)
    X, Y = np.meshgrid(x, y)
    # Iterate through each ellipse and assign values
    for i, (xcen, ang, R1, R2, _, _) in enumerate(data):
        print(xcen, ang, R1, R2)
        xc, yc = xcen
        # Transform coordinates to ellipse frame
        X_rot = (X - xc) * np.cos(-ang) + (Y - yc) * np.sin(-ang)
        Y_rot = -(X - xc) * np.sin(-ang) + (Y - yc) * np.cos(-ang)
        inside_ellipse = (X_rot**2 / R1**2 + Y_rot**2 / R2**2) <= 1
        matrix[inside_ellipse] = nval[i] - 1
    return matrix

def venn_matrix(Ngrid, x1, R1, nval):
    x = np.linspace(-1, 1, Ngrid)
    y = np.linspace(-1, 1, Ngrid)
    X, Y = np.meshgrid(x, y)
    # Initialize matrix
    matrix = np.zeros((Ngrid, Ngrid))
    # Circle centers
    x_left, y_left = -x1, 0
    x_right, y_right = x1, 0
    # Circle equations
    left_circle = (X - x_left) ** 2 + (Y - y_left) ** 2 <= R1 ** 2
    right_circle = (X - x_right) ** 2 + (Y - y_right) ** 2 <= R1 ** 2
    # Intersection area
    intersection = left_circle & right_circle
    matrix[left_circle] = nval[0] - 1  # Left circle
    matrix[right_circle] = nval[1] - 1 # Right circle
    matrix[intersection] = nval[2] - 1 # Intersection area
    return matrix

def bullseye_matrix(Ngrid, nval):
    x = np.linspace(-1, 1, Ngrid)
    y = np.linspace(-1, 1, Ngrid)
    X, Y = np.meshgrid(x, y)
    # Initialize matrix
    matrix = np.zeros((Ngrid, Ngrid))
    # Circle equations
    outer = X ** 2 + Y ** 2 <= 1
    inner = X ** 2 + Y ** 2 <= 0.8 ** 2
    ellipse = X**2/(0.6 ** 2) + (Y+0.0184)**2/(0.4 ** 2) <= 1
    matrix[outer] = nval[0] - 1  
    matrix[inner] = nval[1] - 1
    matrix[ellipse] = nval[2] - 1
    return matrix

def squares_matrix(Ngrid, nval):
    x = np.linspace(-1, 1, Ngrid)
    y = np.linspace(-1, 1, Ngrid)
    X, Y = np.meshgrid(x, y)
    matrix = np.zeros((Ngrid, Ngrid))
    x1, x2, x3 = -0.4, 0.4, 0.6
    y1, y2, y3, y4 = -0.4, -0.2, 0.1, 0.4
    mask_a = (0 <= X) & (X <= x2) & (y2 <= Y) & (Y <= y3)
    mask_b = (x2 <= X) & (X <= x3) & (y2 <= Y) & (Y <= y3)
    mask_c = (x1 <= X) & (X <= x2) & (y1 <= Y) & (Y <= y4)
    matrix[mask_c] = nval[0] - 1
    matrix[mask_b] = nval[1] - 1
    matrix[mask_a] = nval[2] -1
    return matrix
    