import snrTool as snr
import numpy as np 
from matplotlib import pyplot as plt 
from anastruct import SystemElements 
import timeit 


def saveElementModel(ss, saveDir, saveName): 

    folder = snr.createFolder(saveDir, saveName)

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


def getHexagonLatticeCenters(initial_point, lattice_param, x_hexagons, y_hexagons, invertHex=True): 
    
    if invertHex == False: 
    
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

                centerPoints[i, j] = (curr_x, curr_y) 

                # Make sure curr_y is reset 
                if i == x_hexagons-1: 
                    curr_y = start_height 

    else: 

        centerPoints = np.empty( shape=( x_hexagons+1, y_hexagons, 2 ) ) 
        
        # Distance between hexagon center parameters 
        hex_x_distance = 2*lattice_param*np.cos(np.pi/6) # distance from 1 hexagon center to another in x-dir
        hex_y_distance = lattice_param + lattice_param*np.sin(np.pi/6) # distance from 1 hexagon center to another in y-dir 

        # Hexagon grid parameters (from center point to edges) 
        # hex_width = 2*lattice_param 
        # hex_height = 1*np.sqrt(3)*lattice_param

        curr_x = initial_point[0] 
        curr_y = initial_point[1] 

        # initial status of x-shift for each row 
        hex_status = 0 

        for j in range(y_hexagons): 

            if hex_status == 0: 
                curr_x = initial_point[0] 
            elif hex_status == 1: 
                curr_x = initial_point[0] - 0.5*hex_x_distance 
            
            if hex_status == 0: 
                xHexagonCount = x_hexagons 
            elif hex_status == 1: 
                xHexagonCount = x_hexagons + 1 
            for i in range(xHexagonCount): 
                
                if i != 0: 
                    curr_x += hex_x_distance 

                centerPoints[i,j] = (curr_x, curr_y) 
            
            # after processing row of hexagons, switch hex_status variable to shift 
            if hex_status == 0: 
                hex_status = 1 
            elif hex_status == 1: 
                hex_status = 0 
            
            curr_y += hex_y_distance


    return centerPoints 

def hexagonPoints(center, size, roundDecimals=2, invertHex=False): 
    '''
    This function  generates points to throw to anastruct system elements 
    to draw 2D honeycomb stuctures. 

    Args: 
        (center):
            (tuple): Center coordinates of hexagon to be drawn. 
    '''

    points = [] 

    for i in range(6): 

        # Start angel at 0 and increment by 60 each time 
        if invertHex == False: 
            current_angle = 60*(i+1) 
        else: 
            # Start angle at 30 and increment by 60 each time 
            current_angle = 30 + 60*(i) 
        angle_radians = (np.pi/180)*current_angle 

        curr_x = center[0] + size*np.cos(angle_radians) 
        curr_y =  center[1] + size*np.sin(angle_radians) 

        curr_x = np.round(curr_x, roundDecimals) 
        curr_y = np.round(curr_y, roundDecimals) 

        if curr_x == -0.0: 
            curr_x = 0.0 
        if curr_y == -0.0: 
            curr_y = 0.0 

        points.append( ( curr_x, curr_y ) ) 

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

# # --------------------------- Old Node Search Algorithm --------------------------------------- # 

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

def findLowestNodes(ss, startPoint, latticeParam, xHexagons, yHexagons): 
    '''
    Gets lowest node points of the given 2D hexagon crystal lattice. 
    Returns node id's and (x,y) locations for further processing. 

    Args: 

        (ss): 
            (SystemElements): System elements type object. 
        (startPoint): 
            (tuple): (x,y) coordinates of crystal lattice start point. 
        (latticeParam): 
            (int/float): Lattice parameter of hexagons to be drawn. 
        (xHexagons): 
            (int): Number of hexagons in x-direction of lattice. 
        (yHexagons): 
            (int): Number of hexagons in y-direction of lattice. 
    
    Returns: 

        (return1): 
            (list[int, tuple]): List containing node ID of SystemElements object as well as 
            the (x,y) coordinates of each node found. 

    '''

    nodes = [] 
    
    hex_x_distance = 2*latticeParam*np.cos(np.pi/6)
    hex_y_distance = latticeParam + latticeParam*np.sin(np.pi/6)

    curr_x = startPoint[0] 
    curr_y = startPoint[1] - latticeParam 

    for i in range(xHexagons): 

        curr_node_location = ( np.round(curr_x,2), np.round(curr_y,2) ) 
        curr_node_id = ss.find_node_id( curr_node_location )

        try: 
            assert curr_node_id != None
        except AssertionError: 
            print("Assertion Error: No node corresponding to (x,y)={} was found.".format(curr_node_location))
            return None 

        nodes.append( (curr_node_id, curr_node_location) ) 

        curr_x += np.round(hex_x_distance, 2) 

    return nodes  
        
def findHighestNodes(ss, startPoint, latticeParam, xHexagons, yHexagons): 
    '''
    Gets highest node points of the given 2D hexagon crystal lattice. 
    Returns node id's and (x,y) locations for further processing. 

    Node location calculation/iteration is dependent on the number of hexagon
    lattice rows are present. If the number of rows are an odd number, an
    offset of 0.5*hex_x_distance needs to be incremented to the node search
    algorithm to account for the alternating hexagon structure. 

    Args: 

        (ss): 
            (SystemElements): System elements type object. 
        (startPoint): 
            (tuple): (x,y) coordinates of crystal lattice start point. 
        (latticeParam): 
            (int/float): Lattice parameter of hexagons to be drawn. 
        (xHexagons): 
            (int): Number of hexagons in x-direction of lattice. 
        (yHexagons): 
            (int): Number of hexagons in y-direction of lattice. 
    
    Returns: 

        (return1): 
            (list[int, tuple]): List containing node ID of SystemElements object as well as 
            the (x,y) coordinates of each node found. 

    '''

    nodes = [] 

    hex_x_distance = 2*latticeParam*np.cos(np.pi/6)
    hex_y_distance = latticeParam + latticeParam*np.sin(np.pi/6)

    # Start x depending on the number of rows of the lattice
    if yHexagons%2 != 0: 
        curr_x = startPoint[0] 
        nodeCount = xHexagons 
    else: 
        curr_x = startPoint[0] - 0.5*hex_x_distance 
        nodeCount = xHexagons + 1 

    
    # Start y parameter based on number of rows 
    curr_y = startPoint[1] + hex_y_distance*(yHexagons-1) + latticeParam 

    for i in range(nodeCount): 

        curr_node_location = (np.round(curr_x,2), np.round(curr_y,2)) 
        curr_node_id = ss.find_node_id( curr_node_location )

        try: 
            assert curr_node_id != None 
        except AssertionError: 
            print("Assertion Error: No node corresponding to (x,y)={} was found.".format(curr_node_location))

        nodes.append( (curr_node_id, curr_node_location) ) 

        curr_x += np.round(hex_x_distance, 2) 

    return nodes 




def initializeLatticeShell(startPoint, latticeParam, xHexagons, yHexagons, force=10, theta=0, elasticModulus=17000000000,
    moment=5000000, area=0.0005, g=0, invertHex=True, saveDir=None): 
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
    # Initialize finite element material parameters
    EA = elasticModulus*area
    EI = elasticModulus*moment
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
    lattice_coordinates = np.empty(shape= (xHexagons + 1 , yHexagons, 7, 2) )

    # Take not of all elements that are already saved so we don't repeat them 
    # Save initialized system elements points as an array to be referenfced when 
    # creating new system elements to the FEM mesh. 
    old_ss = [] 

    # Initialize hexagon lattice centers: 
    lattice_centers = getHexagonLatticeCenters(startPoint, latticeParam, xHexagons, yHexagons, invertHex=invertHex) 
    print(lattice_centers)

    # Iterate through each hexagon center and fill lattice_coordinates numpy array 
    for j in range(yHexagons): 

        # Every odd numbered row (0 start counting) will have 1 extra hexagon to be drawn in. 
        if j%2 != 0: 
            xHexagonCount = xHexagons + 1 
        else: 
            xHexagonCount = xHexagons 
        
        for i in range(xHexagonCount): 

            # Get hexagon points  
            curr_hexPoints = hexagonPoints(lattice_centers[i,j], latticeParam,invertHex=invertHex) 

            # Initialize start point 
            lattice_coordinates[i, j, 0] = lattice_centers[i,j] 
            lattice_coordinates[i, j, 1:] = curr_hexPoints

            # Iterating through first hexagon being drawn 
            if i == 0: 

                # For even values: 
                if j%2 != 0: 
                    # Finds the point pair list of the right half of the hexagon 
                    point_pair_list = [ (curr_hexPoints[1], curr_hexPoints[0]), 
                                        (curr_hexPoints[0], curr_hexPoints[5]),
                                        (curr_hexPoints[5], curr_hexPoints[4]), ]
                    for k in range(3): 

                        # Check if current element is already added 
                        is_added, index = checkSystemElement( point_pair_list[k][0], point_pair_list[k][1], old_ss ) 

                        # Check if element exists
                        if is_added == False: 
                            ss.add_element( location=[ point_pair_list[k][0], point_pair_list[k][1] ], EA=EA, EI=EI, g=g )
                            old_ss.append( (point_pair_list[k][0], point_pair_list[k][1]) )
                        else: 
                            print("Point Pair {} has already been added at index {} of old_ss list.".format(curr_point_pair, index))

                else: 
                    # Iterate through all point pairs: 
                    
                    for k in range(6): 

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


            elif j%2 != 0 and i == xHexagons: 

                # Uneven row numbers
                if j%2 != 0: 
                    # Finds the point pair list of the right half of the hexagon 
                    point_pair_list = [ (curr_hexPoints[1], curr_hexPoints[2]), 
                                        (curr_hexPoints[2], curr_hexPoints[3]),
                                        (curr_hexPoints[3], curr_hexPoints[4]), ]
                    for k in range(3): 

                        # Check if current element is already added 
                        is_added, index = checkSystemElement( point_pair_list[k][0], point_pair_list[k][1], old_ss ) 

                        # Check if element exists
                        if is_added == False: 
                            ss.add_element( location=[ point_pair_list[k][0], point_pair_list[k][1] ], EA=EA, EI=EI, g=g )
                            old_ss.append( (point_pair_list[k][0], point_pair_list[k][1]) )
                        else: 
                            print("Point Pair {} has already been added at index {} of old_ss list.".format(curr_point_pair, index))

                # Even hexagon values: 
                else: 
                    
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

            # add all system elements and save point pairs 
            else: 

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


    # add supports and flattening elements 
    for i in range(len(lowest_nodes)): 

        ss.add_element( location=[ lowest_nodes[i-1][1], lowest_nodes[i][1] ] )
        old_ss.append( (lowest_nodes[i-1], lowest_nodes[i]) )
        
        # add fixed support on lowest nodes...
        ss.add_support_fixed(node_id=lowest_nodes[i][0])

    # add forces and node locations 
    for i in range(len(highest_nodes)): 

        ss.add_element( location=[ highest_nodes[i-1][1], highest_nodes[i][1] ] )
        old_ss.append( (highest_nodes[i-1], highest_nodes[i]) )

        # add forces at highest nodes...
        # ss.q_load(q=-1, element_id=highest_nodes[i], direction='element') 
        ss.point_load(node_id=highest_nodes[i][0], Fx=0, Fy=force, rotation=theta) 


    # Solve and time it
    solveTimeStart = np.round( timeit.default_timer() , 2 )  
    ss.solve() 
    solveTimeEnd = np.round( timeit.default_timer() , 2 )  
    total_solve_time = solveTimeEnd - solveTimeStart

    print("\n\nSolving beam element model of x-y dimensions, {}-{}, lattice parameter, {},  took {} seconds.\n\n".format(xHexagons, yHexagons, latticeParam, total_solve_time))  


    # Calculate average displacement of nodes
    totalDisplacement = 0 
    displacements = ss.get_node_displacements(node_id=0)
    for i in range(len(displacements)): 
        totalDisplacement += np.sqrt( displacements[i][1]**2 + displacements[i][2]**2 ) 
    averageDisplacements = totalDisplacement/len(displacements) 

    print("\n\nAverage Node Displacement: {}mm".format(averageDisplacements*1000)) 

    if saveDir != None: 
        diagnosticData = [] 
        diagnosticData.append("Hexagon Lattice ({}-x-hexagons,{}-y-hexagons)\n".format(xHexagons, yHexagons))
        diagnosticData.append("Lattice Parameter: {}mm\n".format(latticeParam)) 
        diagnosticData.append("Force: {}N\n".format(force)) 
        diagnosticData.append("Elastic Modulus: {}MPa\n".format(elasticModulus/1000)) 
        diagnosticData.append("Second Moment of Area: {}\n".format(moment)) 
        diagnosticData.append("Area: {}mm\n".format(area*1000))
        diagnosticData.append("Bending Stiffness: {}\n".format(elasticModulus*moment))
        diagnosticData.append("Average Displacement: {}mm\n".format(averageDisplacements)) 
        diagnosticData.append("Mesh Generated in {}s\n".format(total_time))
        diagnosticData.append("Beam Element Model Solved in {}s\n".format(total_solve_time))
    

        folderDir = snr.createFolder( saveDir, "a={}mm A={}mm f={}N E={}MPa".format(latticeParam, area, force, elasticModulus/1000))
        textFile = open("{}\\diagnostics.txt".format(folderDir), 'w') 
        textFile.writelines(diagnosticData) 
        textFile.close() 

        # beam element model is solved in the call of this function. 
        saveElementModel(ss, folderDir, "a={}mm A={}mm f={}N E={}MPa".format(latticeParam, area, force, elasticModulus/1000))          

    return ss 

'''
Genetic Algorithms section: 
1. Initialization of hexagon lattice structure 
2. Fitness measure will equal average displacement values for each node. 
- create function for this 
3. Fix node impact parameters 
- coninuous x-direction testing 
'''
