import snrTool 
import numpy as np 
from matplotlib import pyplot as plt 
from anastruct import SystemElements 
import timeit 

def getHexagonLatticeCenters(initial_point, lattice_param, x_hexagons, y_hexagons): 
    
    # Center points of hexagon lattice to be returned
    centerPoints = np.empty( shape=( x_hexagons, y_hexagons, 2 ) ) 
    
    # Distance between hexagon center parameters 
    hex_x_distance = lattice_param*1.5
    hex_y_distance = np.tan(np.pi/3)*0.5*lattice_param

    # Hexagon grid parameters (from center point to edges) 
    # hex_width = 2*lattice_param 
    # hex_height = 1*np.sqrt(3)*lattice_param

    curr_x = initial_point[0] 
    curr_y = initial_point[1] 

    for j in range(y_hexagons): 

        # Run this code for rows after first row 
        if j != 0: 
            # Initial y-offset 
            curr_y += hex_y_distance*2 
        # Reset curr_x after every row generated
        curr_x = initial_point[0] 
        
        # Initialize start height 
        start_height = curr_y

        for i in range(x_hexagons):

     
            # Initial x-shift 
            curr_x = i*hex_x_distance 

            if i%2 != 0:
                curr_y += hex_y_distance  
        
            else: 
                curr_y = start_height 

            centerPoints[j, i] = (curr_x, curr_y) 

            # Make sure curr_y is reset 
            if i == x_hexagons-1: 
                curr_y = start_height 

    return centerPoints 

def hexagonPoints(center, size, horizontal=True, roundDecimals=2): 
    '''
    This function  generates points to throw to anastruct system elements 
    to draw 2D honeycomb stuctures. 

    Args: 
        (center):
            (tuple): Center coordinates of hexagon to be drawn. 
    '''

    if horizontal == False: 
        points = [] 

        for i in range(6): 
            current_angle = 60*(i+1) - 30
            angle_radians = (np.pi/180)*current_angle 

            curr_x = center[0] + size*np.cos(angle_radians) 
            curr_y = center[1] + size*np.sin(angle_radians) 

            points.append((curr_x, curr_y)) 

    else: 

        points = [] 

        for i in range(6): 
            current_angle = 60*(i+1) 
            angle_radians = (np.pi/180)*current_angle 

            curr_x = center[0] + size*np.cos(angle_radians) 
            curr_y =  center[1] + size*np.sin(angle_radians) 

            points.append( (np.round(curr_x, roundDecimals), np.round(curr_y, roundDecimals) ) ) 

    return points 

def checkSystemElement(point1, point2, pointPairList): 
    '''
    This function searches through a point list to see if a pair 
    of points is present already in the pointsList. To be used so that
    an anaStruct object of class SystemElements() does not repeat an element. 
    This fucks up the total FE simultion. 

    Args:

        (point1): 
            (tuple): First point in the element. 
        (point2):
            (tuple): Second point in the element. 
        (pointPairList): 
            (list): List of points to be searched. 

    Returns: 
        (return1): 
            (boool): Whether or not the point pair has been found 

        (return2): 
            (int/None): If the point pair is found within the pointList, the index
            at which the point pair is found will be returned. Otherwise if the point pair
            is not found, a None type will be returned. 
    '''
    try: 
        index = pointPairList.index( (point1, point2) )
    except ValueError: 
        
        # If not found, try different orientation of points 
        try: 
            index = pointPairList.index( (point2, point1) ) 
        except ValueError: 

            # Pair of points consisting of point1 and point2 do not
            # exist in the given pointPairList[]. 
            return False, None  
    
    # if it doesn't return false in the above code, return the index
    return True, index 


'''
I messed up this function...wrote it the lazy way (hard coding xddd) 
Check findLowestNodes() function below it. 
'''
# def findLowestNodes(xHexagons, yHexagons): 
#     '''
#     This function finds the node-id of the lowest nodes in the x-y cartesian 
#     plane of which matplotlib has drawn the elements. 
#     '''
    
#     if (xHexagons%2 != 0): 
#         num_nodes = (xHexagons+1)//2 
#     else:
#         num_nodes = xHexagons//2

#     node_ids = [] 

#     for i in range(num_nodes):

#         node_ids.append(i*10 + 5)
#         node_ids.append(i*10 + 6) 

#     return node_ids 

# --------------------------- Old Node Search Algorithm --------------------------------------- # 

# def findLowestNodes(ss, startPoint, latticeParam, xHexagons, yHexagons): 

#     # hex_x_distance = latticeParam*1.5
#     hex_y_distance = np.tan(np.pi/3)*0.5*latticeParam

#     # x distance when hexagon changes height 
#     hex_x_short = np.cos(np.pi/3)*latticeParam 

#     # Initialize x, y positions of bottom left most node. 
#     curr_y = np.round( startPoint[1] - hex_y_distance , 2 )
#     curr_x = np.round( startPoint[0] - 0.5*latticeParam , 2 ) 

#     node_ids = [] 

#     count = 0 
#     node_id = -1 

#     for i in range(xHexagons*2): 

#         if i%2 == 0: 
#             # curr_y = np.round( startPoint[1] + (yHexagons*2 - 1)*hex_y_distance , 2 )
#             x_step = latticeParam 
#         else: 
#             # curr_y = np.round( startPoint[1] + (yHexagons*2)*hex_y_distance , 2 )
#             x_step = hex_x_short 

#         if node_id == -1: 
#             curr_y = np.round( startPoint[1] - hex_y_distance , 2 )
#         else: 
#             curr_y = startPoint[1]

#         curr_node_id = ss.find_node_id( (curr_x, curr_y) )

#         node_ids.append(curr_node_id)

#         # increment curr_x with x_step 
#         curr_x += x_step 

#         count += 1 

#         if count == 2: 
#             node_id = node_id*(-1) 
#             count = 0 

#     return node_ids 

# def findHighestNodes(ss, startPoint, latticeParam, xHexagons, yHexagons): 

#     # hex_x_distance = latticeParam*1.5
#     hex_y_distance = np.tan(np.pi/3)*0.5*latticeParam

#     # x distance when hexagon changes height 
#     hex_x_short = np.cos(np.pi/3)*latticeParam 

#     # Initialize x, y positions of top left most node. 
#     curr_y = np.round( startPoint[1] + (yHexagons*2 - 1)*hex_y_distance , 2 )
#     curr_x = np.round( startPoint[0] - 0.5*latticeParam , 2 ) 

#     node_ids = [] 

#     count = 0 
#     node_id = -1 

#     for i in range(xHexagons*2): 

#         if i%2 == 0: 
#             # curr_y = np.round( startPoint[1] + (yHexagons*2 - 1)*hex_y_distance , 2 )
#             x_step = latticeParam 
#         else: 
#             # curr_y = np.round( startPoint[1] + (yHexagons*2)*hex_y_distance , 2 )
#             x_step = hex_x_short 

#         if node_id == -1: 
#             curr_y = np.round( startPoint[1] + (yHexagons*2 - 1)*hex_y_distance , 2 )
#         else: 
#             curr_y = np.round( startPoint[1] + (yHexagons*2)*hex_y_distance , 2 )

#         curr_node_id = ss.find_node_id( (curr_x, curr_y) )

#         node_ids.append(curr_node_id)

#         # increment curr_x with x_step 
#         curr_x += x_step 

#         count += 1 

#         if count == 2: 
#             node_id = node_id*(-1) 
#             count = 0 

#     return node_ids 

# --------------------------- Old Node Search Algorithm --------------------------------------- # 


def initializeLatticeShell(startPoint, latticeParam, xHexagons, yHexagons, force, theta=0, EA=15000
    , EI=5000, g=0, flushBottomTop=True): 
    '''
    This function takes in initial parameters for a hexagon lattice, creates
    a mesh structure using the anaStruct library object SystemElements(). First, 
    the center points of hexagons in the lattice will be generated, the hexagon 
    vertices will be saved in a custom NumPy array to save the center points and 
    al 6 vertices. 

    The numpy array will hold 7 tuples( 1 tuple for the center-point, 6 tuples for
    coordinates of each hexagon vertex. )

    Args: 

        (startPoint): 
            (tuple): Tuple of dtype=np.float64 or dtype=np.float32 housing the 
            original start point (bottom left) of the hexagon lattice to be constructed.
        (latticeParam): 
            (float): The length of each 'sub triangle' length starting from the center
            of the hexagon to each vertex. 
        (xHexagons): 
            (int): Number of hexagons to be constructed in the x-direction in the lattice.
        (yHexagons): 
            (int): Number of hexagons to be constructed in the y-direction in the lattice.
    
    Returns: 

        (lattice): 
            (SystemElements() Object): The system elements object containing all trusses
            defined and fully constrained. 
        (lattice_coordinatews): 
            (numpy array): Custom numpy array of shape=( xHexagons, yHexagons, 7, 2 ). 
            Contains 7 coordinates per hexagon: Center coords, and vertex coordinates.
    '''

    # Initialize timer 
    start = timeit.default_timer() 

    # Initialize anaStruct.py SystemElements() object...
    # ss = SystemElements(figsize=(12,8), EA=15000.0, EI=5e3) 
    ss = SystemElements(EA=EA, EI=EI) 
    '''
    EA - Standard axial stiffness of element 
    EI - Standard bending stiffness of elements 
    figsize - Matplotlibs standard figure size (inches) 
    '''

    # Initialize custom numpy array...
    lattice_coordinates = np.empty(shape= (xHexagons, yHexagons, 7, 2) )

    # Take not of all elements that are already saved so we don't repeat them 
    # Save initialized system elements points as an array to be referenfced when 
    # creating new system elements to the FEM mesh. 
    old_ss = [] 

    # Initialize hexagon lattice centers: 
    lattice_centers = getHexagonLatticeCenters(startPoint, latticeParam, xHexagons, yHexagons) 
    print(lattice_centers)

    # Iterate through each hexagon center and fill lattice_coordinates numpy array 
    for i in range(yHexagons): 
        
        for j in range(xHexagons): 

            # Get hexagon points 
            curr_hexPoints = hexagonPoints(lattice_centers[i,j], latticeParam) 

            # Initialize start point 
            lattice_coordinates[i, j, 0] = lattice_centers[i,j] 
            lattice_coordinates[i, j, 1:] = curr_hexPoints

            # add all system elements and save point pairs 
            for k in range(6): 

                # check if current element is already added...if not, add it! 
                curr_point_pair = ( curr_hexPoints[k-1], curr_hexPoints[k] )
                is_added, index = checkSystemElement( curr_point_pair[0], curr_point_pair[1], old_ss )

                # If element doesn't already exist, add it to the mesh 
                if is_added == False: 

                    ss.add_element( location=[ curr_hexPoints[k-1], curr_hexPoints[k] ], EA=EA, EI=EI, g=g )
                    old_ss.append( (curr_hexPoints[k-1], curr_hexPoints[k]) )

                else: 
                    
                    # save elements added 
                    # old_ss.append( curr_hexPoints[k-1], curr_hexPoints[k] )

                    # testing line --------------------
                    print("Point Pair {} has already been added at index {} of old_ss list.".format(curr_point_pair, index))
                    # end testing line ------------------


    end = timeit.default_timer() 

    total_time = np.round( (end-start), 2 ) 

    print("\n\nGenerating hexagon mesh of x-y dimensions, {}-{}, lattice parameter, {},  took {} seconds.\n\n".format(xHexagons, yHexagons, latticeParam, total_time))  

    lowest_nodes = findLowestNodes(ss, startPoint, latticeParam, xHexagons, yHexagons) 
    highest_nodes = findHighestNodes(ss, startPoint, latticeParam, xHexagons, yHexagons) 

    # add supports here... 
    for i in range(len(lowest_nodes)): 
        
        # add fixed support on lowest nodes...
        ss.add_support_fixed(node_id=lowest_nodes[i])

    # add forces here... 
    for i in range(len(highest_nodes)): 

        # add forces at highest nodes...
        # ss.q_load(q=-1, element_id=highest_nodes[i], direction='element') 
        ss.point_load(node_id=highest_nodes[i], Fx=0, Fy=force, rotation=theta) 


    return ss 

def saveElementModel(ss, saveDir, saveName): 

    folder = snrTool.createFolder(saveDir, saveName)

    fig = ss.show_axial_force(show=False)
    plt.title("axial force")
    plt.savefig("{}\\{}.png".format(folder,"axial force")) 

    fig = ss.show_bending_moment(show=False)
    plt.title("bending moment")
    plt.savefig("{}\\{}.png".format(folder, "bending moment")) 

    fig = ss.show_displacement(show=False) 
    plt.title("displacement")
    plt.savefig("{}\\{}.png".format(folder, "displacement")) 

    fig = ss.show_reaction_force(show=False) 
    plt.title("reaction force")
    plt.savefig("{}\\{}.png".format(folder, "reaction force")) 

    fig = ss.show_shear_force(show=False)
    plt.title("shear force")
    plt.savefig("{}\\{}.png".format(folder, "shear force")) 

    fig = ss.show_structure(show=False) 
    plt.title("structure")
    plt.savefig("{}\\{}.png".format(folder, "structure")) 


'''
Genetic Algorithms section: 
1. Initialization of hexagon lattice structure 
2. Fitness measure will equal average displacement values for each node. 
- create function for this 
3. Fix node impact parameters 
- coninuous x-direction testing 
'''