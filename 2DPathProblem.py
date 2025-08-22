Problem Description: 

Write a Python function that determines whether a path exists between two points that only crosses black pixels. Define a path as a sequence of pixels such that each pixel in the path is adjacent to the next pixel in the path. 

Solution Description: 

This problem is solved as an X Y axis traversal problem. 

User defined values are 
  * the 2D Environment size: NxM 
    example:[4,4] is a square N=4, M=4 matrix of random 0,1 values 
  * start_position within array = [x1,y1] 
  * end_position within array = [x2,y2]
  * flag_value = 1

Given user defined inputs of the 2D Environment size and start and end positions, the program will create a 2D Enviroment of the user defined values N x M with random 0,1 values. The axis distances across the X and Y axis are determined in each axis direction. 

Path Determination Approach: First, a function, axis_steps, takes the difference in the axis movements: Y2-Y1 and X2-X1 to determine 
axis distance ranges. Using the distances, a dictionary of step movements are created for both the X and Y axis. 
A function, build_xy_path_dic, is used to iterate through the two axis step dictionaries using a larger one if exists as at the bound and populates 
different path potentials from A to B. This is done by a permutation approach of adding X and Y and then Y then X to each previous step. 

There is also a function, produce_xmax_ymax_path_ymax_max_path, to produce the maximum path first along the X axis 
and then along the Y axis and again, first along the Y axis and then along the X axis. This will be used in combination 
with the XY combination paths for a complete path set. 

All possible paths produced are combined to generate the list of the possible paths of index values. Finally, out of the list of possibilies
of index value paths, it is determined if there exist one of only the flag value, i.e a path of all 1's.

import numpy as np

# User Defined values of 2D Enviroment size, the start and end positions and flag_value 
N = 7
M = 7
start_position = [2,2]
end_position = [5,5]
flag_val = 1


# Build out matrix
array = np.random.choice([1,1],[N,M])


'''
Funcion Description: Produces a step dictionary in both X and Y directions
Input: start_position, end_position 
Output: x_dic_steps, y_dic_steps 
'''

def produce_xy_step_dic(start_position, end_position):
    # Determine difference in X and Y steps needed 
    start_position_x = start_position[0]
    start_position_y = start_position[1] 

    end_position_x = end_position[0]
    end_position_y = end_position[1]

    
    # Number of Steps: in the X axis and Y axis directions  
    n_x_steps = end_position_x - start_position_x 
    n_y_steps = end_position_y - start_position_y 
    
    # Build out dictionary of X steps and Y steps 
    x_dic_steps = {}
    y_dic_steps = {}
    xy_dics = [] 
    for i in range(0,n_x_steps):
        i=i+1
        pos_i_x = [start_position_x + i , start_position_y]
        x_dic_steps[i] = pos_i_x
    
    for i in range(0,n_y_steps):
        i=i+1
        pos_i_y = [start_position_x, start_position_y+i]
        y_dic_steps[i] = pos_i_y

    xy_dics = [x_dic_steps, y_dic_steps]
    
    return xy_dics


'''
Function Description: Intakes a dictionary and outputs an array of steps in order to be applied to 
    produce:
        X_MAX,Y_MAX paths
        Y_MAX,X_MAX paths 

Input: Dictionary of steps
Output: Iterated list of steps along an axis direction  
'''

def axis_steps(dic):

    size_dic = len(dic.keys())
    i_pos_path = [] 
    
    for i in range(0,size_dic):
        i=i+1
        i_position = dic[i]
        i_pos_path.append(i_position)
        
    return i_pos_path



'''
Function Description: Intakes both X and Y dictionaries to produce two paths:
    XMAX_YMAX = The path formed when traversing first the maximum path value along only the X axis first followed by Y traversal 
    YMAX_XMAX = The path formed when traversing first the maximum path value along only the Y axis first followed by X traversal 

    Both will be used in final path possible combinations and appened to paths when switching between X and Y axis moves 
    Note: Total paths should be 2 + N_steps*2 = max path in x then y + max path in y then x + N_steps*2 where 2 is the choice selection of X or Y

Input: X step dictionary, Y step dictionary, axis_step function
Output: List containing two paths: [XMAX_YMAX, YMAX_XMAX]
'''

def produce_xmax_ymax_path_ymax_max_path(dic_x, dic_y, axis_steps):

    x_max = []
    y_max = [] 
    
    y_path = axis_steps(dic_y)
    x_path = axis_steps(dic_x)  
    
    
    for i in x_path:
        x_max.append(i)
    
    for j in y_path:
        y_max.append(j)
    
    xy_max = [] 
    yx_max = [] 
    
    for i in x_max:
        xy_max.append(i)
        
    for j in y_max:    
        xy_max.append(j)
        yx_max.append(j)
        
    for k in x_max:
        yx_max.append(k)

    xy_max_yx_max_paths = [
                            xy_max,
                            yx_max
                            ] 
    return xy_max_yx_max_paths


'''
 Function Description: Building out the XY combination path dictionary
          Logic takes the larger of the two dicionaries to determine key path set up
          Iterate through the larger dictionary up to the size of the smaller
          Then proceed forward only with the larger dictionary values adding on to the built out dictionary in the preceeding steps

 
 Input: X Step Dictionary,  Y Step Dictionary 
 Output: Dictionary of the different XY path combinations from start_position to end_position
         This will need to be combined with the two paths that comprise going only along the maximum direction of X then Y and Y then X
         

Example:

For each dictionary key, store the combination set of appending 
                            an xy_step to all previous lists 
                            an yx_step to all previous lists

path_dic_1_step = {[1]:  [[0,1],[1,0]],
                         [[1,0],[0,1]]
                }


            
path_dic_2_Steps           {[1]: [[0,1],[1,0]],
                                 [[1,0],[0,1]],
                                 
                            [2]: [[0,1],[1,0],[0,2],[2,0]],
                                 [[1,0],[0,1],[0,2],[2,0]],
                                 [[0,1],[1,0],[2,0],[0,2]],
                                 [[1,0],[0,1],[2,0],[0,2]]
                    
                            }

path_dic_2_Steps           {[1]: [[0,1],[1,0]],
                                 [[1,0],[0,1]],
                                 
                            [2]: [[0,1],[1,0],[0,2],[2,0]],
                                 [[1,0],[0,1],[0,2],[2,0]],
                                 [[0,1],[1,0],[2,0],[0,2]],
                                 [[1,0],[0,1],[2,0],[0,2]],

                            # all of step two values appending xy_i step and all of step values appending yx_i step
                            [3]:
                             
                    
                            }

            
'''


def build_xy_path_dic(y_dic_steps, x_dic_steps):
    
    size_steps_y = len(y_dic_steps.keys())
    size_steps_x = len(x_dic_steps.keys())
    
    paths = {}
    path_dic = {}
    
    
    if size_steps_y > 0 and size_steps_x > 0:
        if size_steps_y >= size_steps_x:  
            for i in range(0, size_steps_y):
                i=i+1
                
                if i <= size_steps_x:
                    a_step_i = y_dic_steps[i]
                    b_step_i = x_dic_steps[i]
                    
                    ab_i = [a_step_i] + [b_step_i]
                    ba_i = [b_step_i] + [a_step_i]
                    
                    if len(path_dic.keys()) == 0:
                        path_dic[i] = [ab_i,ba_i]
                                      
                    else:
                        path_lists = path_dic[i-1] 
                    
                        # intialize empty new path list to combine lists and store in dic[i+1] key 
                        new_list = [] 
                        for j in path_lists:
                            new_list.append(j + ab_i)
                        for j in path_lists:
                            new_list.append(j + ba_i)
                    
                        path_dic[i] = new_list
                
                else:
                    # Then we have  j > i so only dealing with larger dic just look at  
                    y_j = y_dic_steps[i]
                    
                    path_lists = path_dic[i-1]
                    new_list = [] 
                    for path in path_lists:
                        new_list.append(path + [y_j])
                    
                    path_dic[i] = new_list
            
    if size_steps_y == 0 and size_steps_x > 0:
        for i in range(size_steps_x):
            i=i+1
            a_step_i = 0
            b_step_i = x_dic_steps[i]
                    
            ab_i = [a_step_i] + [b_step_i]
            ba_i = [b_step_i] + [a_step_i]
                    
            if len(path_dic.keys()) == 0:
                path_dic[i] = [ab_i,ba_i]
                                      
            else:
                path_lists = path_dic[i-1] 
                    
            # intialize empty new path list to combine lists and store in dic[i+1] key 
                new_list = [] 
                for j in path_lists:
                    new_list.append(j + ab_i)
                for j in path_lists:
                    new_list.append(j + ba_i)
                    
                path_dic[i] = new_list

    if size_steps_y > 0 and size_steps_x == 0:
        for i in range(size_steps_y):
            i=i+1
            a_step_i = y_dic_steps[i]
            b_step_i = 0
                    
            ab_i = [a_step_i] + [b_step_i]
            ba_i = [b_step_i] + [a_step_i]
                    
            if len(path_dic.keys()) == 0:
                path_dic[i] = [ab_i,ba_i]
                                      
            else:
                path_lists = path_dic[i-1] 
                    
            # intialize empty new path list to combine lists and store in dic[i+1] key 
                new_list = [] 
                for j in path_lists:
                        new_list.append(j + ab_i)
                for j in path_lists:
                        new_list.append(j + ba_i)
                    
                path_dic[i] = new_list
                
               
           
    return path_dic

'''
Function Description: Produces a path of index postition values given a start and end
Input: Start_position, end_position, functinos to produce paths
Output: List of Lists of Path values for the index position
'''

def produce_all_index_value_paths(start_position, end_position, produce_xy_step_dic, produce_xmax_ymax_path_ymax_max_path, build_xy_path_dic):
        
        # For first start,end produce xy steps and path dictionary 
        # we want a path of all 1's 

        
        xy_dics = produce_xy_step_dic(start_position, end_position)
        x_dic_steps = xy_dics[0]
        y_dic_steps = xy_dics[1]
        
        xy_max_yx_max_paths = produce_xmax_ymax_path_ymax_max_path(x_dic_steps, y_dic_steps, axis_steps)
        
        n_x_steps = len(x_dic_steps.keys())
        n_y_steps = len(y_dic_steps.keys())
        
        if n_x_steps >= n_y_steps:
            xy_path_dic = build_xy_path_dic(x_dic_steps , y_dic_steps)
            
        if n_y_steps >  n_x_steps:
            xy_path_dic = build_xy_path_dic(y_dic_steps, x_dic_steps)
        
        
        # Last key is the full list of combo paths
        size_xy_path_dic = len(xy_path_dic.keys())
        
        paths = xy_path_dic[size_xy_path_dic]

        index_value_paths = []

        '''
        N paths of XY -> YX combintions formations
        '''
        for i in range(0,len(paths)):
            path_i = paths[i]
            path_i.insert(0,start_position) 
            path_i.append(end_position)
            index_value_paths.append(path_i) 
        

        '''
        Two paths of maximum Y movement then maximum X movement 
        '''
        xy_max_path = xy_max_yx_max_paths[0]
        yx_max_path = xy_max_yx_max_paths[1]
        
        xy_max_path.insert(0,start_position)
        xy_max_path.append(end_position)
        
        yx_max_path.insert(0,start_position)
        yx_max_path.append(end_position)

        index_value_paths.append(xy_max_path)
        index_value_paths.append(yx_max_path)



        return index_value_paths


'''
Function Description: Intakes a list of indexes and numpy array to produce a list of values per index within the array
                      This will be used to determine if there exist a possible path from start to end of only the flag value 

Input: Path_list, 2D Array Environment 
Output: List of index values 


Function Update: If we want to produce path value that does not intersect another path and it's associated array values 
                then we want to store not just associated array values of the position index but the position index itself as well

                possible: dictionary  [pos_index] = [index_val]
                OR  
                          dictionary  key =  [[index_val],[pos_index]]
'''
def produce_index_path_values(path_list, array):
    
    all_path_dic = {} 
    all_paths = []
    n_index_path_possibilities = len(path_list)
    
    for i in range(0,n_index_path_possibilities):
         # will be dic key so we can get this value if all unique and then see if all associated flag_value: i.e 1's
        index_position_path_list = [] 
        array_path_values = [] 
        index_position_path = path_list[i]
        index_position_path_list.append(index_position_path)
        
        for index_position_val in index_position_path:
            index_0_val = index_position_val[0]
            index_1_val = index_position_val[1]
            array_index_value = array[index_0_val, index_1_val] 
            array_path_values.append(array_index_value)
        
        all_paths.append(array_path_values)

        all_path_dic[str(index_position_path_list)] = array_path_values
    
    return [all_paths, all_path_dic]





'''
Function Description: Produces the pixel position index value based on the X and Y step direction movements
Input: List of lists of the index_paths as a result of X step and Y step movements 
Output: The index position value updated at each step along the path resulting in the path of index possibilites based on
        different permutations of axis movements either along the x and y axis from start point to end 
'''

def produce_index_position_list(index_path_list): 
    index_positions_list = [] 
    
    for path_list in index_path_list:
        # Everytime first index value of each list within the LIst of lists is the start point
        start_val = path_list[0]
        
        start_val_x = start_val[0]
        start_val_y = start_val[1]
              
        size_path = len(path_list)
    
        index_pos_x = start_val_x
        index_pos_y = start_val_y
    
        p_val = [] 
        p_val.append([start_val_x,start_val_y])
        for i in range(1, size_path-1):
                
                path_i_val = path_list[i]
                path_i_val_x = path_i_val[0]
                path_i_val_y = path_i_val[1]
            
                # Observe difference between the _x and _y each step 
                # if delt_x  == 1 then add 1 to x 
                # if delt_y  == 1 then add 1 to y
               
                delt_x = path_i_val_x - index_pos_x
                delt_y = path_i_val_y - index_pos_y
    
                # update start_x start_y 
                if delt_x  == 1:
                    index_pos_x = index_pos_x + delt_x
                if delt_y == 1:
                    index_pos_y = index_pos_y + delt_y
    
                p_val.append([index_pos_x, index_pos_y])
    
            
        index_positions_list.append(p_val) 

        
    return index_positions_list 


'''
Function Description: Intakes two lists and compares; Finds a path where no elemet is within the other list
                      Returns the first matched comparison


Approach: First item from the first list is going to be compared against every item in the second list. However, 
          we really only will compare up to the first found list without an identical value. 

          A = [[a1, a2, a3], [a4,a5,a6]]
          B = [[b1,b2, b3] , [b4,b5,b6]]

          -a1 vs b1
          -a1 vs b2 
          -a1 vs b3
          -a2 vs b1
          -a2 vs b2
          -a2 vs b3
          -a3 vs b1
          -a3 vs b2
          -a3 vs b3

          logic: 
                    if a1 = b1: then stop, move on to next i +1 
                    if a1 != b1: check if 
                 
          
Input: Two List of lists respresenting the index positions: index_position_list_1, index_position_list_2
Output: First list from first and first list from second where there is a unique set of values per list 
'''

def produce_unique_list_pair(lists_1, lists_2):
    
    # For this list 
    # we would want to keep it if we find that comparing each item in it to every item in the another each item is not the same
            
        # Check every item of list_i against every item of list_j
        # if every item list_i != list_j 
        # put the list_i and list_j in the dictionary 
        #    : how:   
        #    :    ap.1  
        #           append each item in list_i to new list 
        #           if size_new_list_i == size_list_i and size_new_list_j == size_list_j
        #           then new appended lists of the unequeal values are the same size as the original list of values checked 
        #           store both new_list_i and new_list_j in the diciontary 


        #    :    ap. 2
        #            keep track of a boolean array indicating if match or notmatch of item_i vs. item_j
        #            if bool_array is all false then return the two lists checked 
        #            if item_j == item_i
        #                set boolean to 1 and move on to next list to check 
        #                i = i+1
        #            if item_j != item_i
        #                bool_val = 0 
        #                if bool_val = 0 and  

    
    # possibilities: 
    #        everything the same so no unmatched set
    #        there exist a list in one that is not the same as a list in the other 
    # input1: list_i = [[a1, a2,a3], [a4,a5,a6]]
    # input2: list_j = [[b1,b2,b3],[b4,b5,b6]
    
    size_lists_1 = len(lists_1)
    size_lists_2 = len(lists_2)
        
        
    dic = {}
        
    # return list_i if it does not match anything in list_j; return list_i, list_j 
    # if item in list_i == item in list_j: move on i = i + 1 
    
    for i in range(size_lists_1):
        
        # empty dictionary so find a qualifiying list to input into dic
        # qualifying list: all items in list1 are NOT in list2
        
        if len(dic.keys()) == 0:
    
            list_i = lists_1[i]
            size_list_i = len(list_i)
            all_val_dic = {} 
            for n in range(size_list_i):
                list_i_item_n = list_i[n]
                
                # Compare to second list set
                for j in range(size_lists_2):
                    
                    # given a list, check if each item in the list, list_i_item_n, matches each item in another list
                    list_j = lists_2[j]
                    size_list_j = len(list_j)
                    
                    # Keep track of all matched values in a list
                    # if all False then store in dictionary 
    
                    val_list =[] 
                    for m in range(size_list_j):
                        list_j_item_m = list_j[m]
                        val = list_i_item_n == list_j_item_m
                        val_list.append(val)
                        all_val_dic[str(list_i_item_n)] = val_list
                                
                
            if all(set(v) == {False} for v in all_val_dic.values()):
                dic['qualifying_lists'] = [list_i , list_j]
    
        else:
            return dic
                


    

'''
Given start_position, end_position, produce both dictionaries of steps in the X direction,Y direction
     resulted by going along: 
                 maximum x direction then y,
                 maximum y direction then x,
                 alterate across first x,y then y,x  or first y,x then x,y
'''

xy_dics = produce_xy_step_dic(start_position, end_position)
x_dic_steps = xy_dics[0]
y_dic_steps = xy_dics[1]


xy_max_yx_max_paths = produce_xmax_ymax_path_ymax_max_path(x_dic_steps, y_dic_steps, axis_steps)


n_x_steps = len(x_dic_steps.keys())
n_y_steps = len(y_dic_steps.keys())

'''
Build out path dictionary of different step combintation in x,y direction  
'''
if n_x_steps >= n_y_steps:
    xy_path_dic = build_xy_path_dic(x_dic_steps , y_dic_steps)
    
if n_y_steps >  n_x_steps:
    xy_path_dic = build_xy_path_dic(y_dic_steps, x_dic_steps)


# Last key is the full list of combo paths
size_xy_path_dic = len(xy_path_dic.keys())
paths = xy_path_dic[size_xy_path_dic]

'''
Given the possible path combinations add the start and end position to get the complete index value path 

Store all path possiblities in a list index_value_paths to produce the index position
'''
index_value_paths = []

for i in range(0,len(paths)):
    path_i = paths[i]
    path_i.insert(0,start_position) 
    path_i.append(end_position)
    index_value_paths.append(path_i) 


xy_max_path = xy_max_yx_max_paths[0]
yx_max_path = xy_max_yx_max_paths[1]

xy_max_path.insert(0,start_position)
xy_max_path.append(end_position)

yx_max_path.insert(0,start_position)
yx_max_path.append(end_position)

index_value_paths.append(xy_max_path)
index_value_paths.append(yx_max_path)

'''
Produce list of the index_position as we move along the the path of x y steps stored in the index value paths produced from x or y moves
'''
index_position_path_list = produce_index_position_list(index_value_paths)

'''
Given the array and the index positions, produce an list of the indexed values per position
This will produce the path whether it is all flag_value or not
'''
position_index_value = produce_index_path_values(index_position_path_list , array)
position_index_value_list = position_index_value[0]



''' 
Check to see if all values within a path are the flag_value
'''
path_all_flag_value = [] 

for index_position_values_path in position_index_value_list:
    if all(v == flag_val for v in index_position_values_path):
        path_all_flag_value.append(index_position_values_path)

total_paths = len(index_position_path_list)
total_paths_only_all_flag_value = len(path_all_flag_value)

print("Total paths:", total_paths)
print("Total Paths Only Flag Value:", total_paths_only_all_flag_value, "Flag_value:", flag_val) 


if total_paths_only_all_flag_value > 0:
    print("True")
    print("Total possible paths all flag_value:",total_paths_only_all_flag_value)
    print("Total Paths Possible:",total_paths)
else:
    print("False", total_paths_only_all_flag_value)
    print("Total Paths POssible:", total_paths) 



# Check of the paths of index position values 
for index_values_path in index_position_path_list:
    print(index_values_path)


# Check of 2D Envirnment 
array

'''
=====================================================================================

Additional Goal: Given a set of two startandend points: two start and two end points
    
    a. Determine if exist a path between both that are all flag_value
    b. determine if exist a path in one both not the other i.e there is 
=====================================================================================

First we can check if both start and end value pairs:
    1. are both the flag values; otherwise stop and return False because there will only exist any path between start and end all flag value 
    if both end and start are also the flag value 
    2. start and end pairs can't be same 

    If both checks are passed then we can proceed forward with computing paths possibilites from start-to-end

    Produce dictionary for both paths 

    Compare across

    VS. Produce dictionary for first path and only produce dicionary for second point if not path in the first dicionary 
'''



start_position_2 = [0,0]
end_position_2 = [3,3]
start_position_1 = [4,4]
end_position_1 = [6,5]
flag_value = 1 

# Given two pairs of points: 
# produce two paths one for each pair:all black and no pixel in both (each has unique index positions)  

# First check to see if both array value of all points are black or flag value 

start_pos_1_index_value = array[start_position_1[0], start_position_1[1]] 
end_pos_1_index_value = array[end_position_1[0], end_position_1[1]] 

start_pos_2_index_value = array[start_position_2[0], start_position_2[1]] 
end_pos_2_index_value = array[end_position_2[0], end_position_2[1]] 


# First check: No path within another so end points can't be equal
#            : Return False if any point is equal; unable to produce unique paths 

if (start_position_2 != start_position_1 and start_position_2 != end_position_1 and end_position_1 != start_position_1 and end_position_1 != end_position_2):
    
    # Second check:  End points are flag value 
    #            :  Return false if any end point is not flag_value
    
    start_pos_1_index_value = array[start_position_1[0], start_position_1[1]] 
    end_pos_1_index_value = array[end_position_1[0], end_position_1[1]] 
    
    start_pos_2_index_value = array[start_position_2[0], start_position_2[1]] 
    end_pos_2_index_value = array[end_position_2[0], end_position_2[1]] 

    
    # If either NOT flag value then there will NOT exist an path between both of all VALUE value
    if (start_pos_1_index_value == flag_value and end_pos_1_index_value == flag_value and start_pos_2_index_value == flag_value and end_pos_2_index_value == flag_value):
        # produce first dic of flag_values start_end position 

        '''
         #Produce list of the index_position as we move along the the path of x y steps stored in the index value paths produced from x or y moves
         
         index_position_path_list = produce_index_position_list(index_value_paths)
         
         
         #Given the array and the index positions, produce an list of the indexed values per position
         #This will produce the path whether it is all flag_value or not
         
         position_index_value_list = produce_index_path_values(index_position_path_list , array)
        
        '''
        path_step_dic_1 = produce_all_index_value_paths(start_position_1, end_position_1, produce_xy_step_dic, produce_xmax_ymax_path_ymax_max_path, build_xy_path_dic)
        path_step_dic_2 = produce_all_index_value_paths(start_position_2, end_position_2, produce_xy_step_dic, produce_xmax_ymax_path_ymax_max_path, build_xy_path_dic)
        
        index_positions_list_1 = produce_index_position_list(path_step_dic_1)
        index_positions_list_2 = produce_index_position_list(path_step_dic_2)


        position_index_value_1 = produce_index_path_values(index_positions_list_1, array)
        position_index_value_2 = produce_index_path_values(index_positions_list_2, array)
        
        position_index_value_1_list = position_index_value_1[0]
        position_index_value_1_dic = position_index_value_1[1]
        
        position_index_value_2_list = position_index_value_2[0]
        position_index_value_2_dic = position_index_value_2[1]


        unique_lists = produce_unique_list_pair(index_positions_list_1, index_positions_list_2)



        qualifying_path_1_values = position_index_value_1_dic[str([unique_lists['qualifying_lists'][0]])]
        qualifying_path_2_values = position_index_value_2_dic[str([unique_lists['qualifying_lists'][1]])]
        
        
        if all(v == flag_val for v in qualifying_path_1_values) and all(v == flag_val for v in qualifying_path_1_values):
            print("True")
            print("dictionary of unique qualifying lists:")
            print(unique_lists['qualifying_lists'])    
        else:
             print("False", (size_index_path_1, size_index_path_2 ))
         
    
         
    else:
             print("False: Path ends does not contain flag_value")
             print("start_pos_1_index_value:",start_pos_1_index_value)
             print("end_pos_1_index_value:",end_pos_1_index_value)
             print("start_pos_2_index_value:",start_pos_2_index_value)
             print("end_pos_2_index_value:",end_pos_2_index_value)
                   
         
        
else:
    print("End points same: No path exists where pixels are unique to both since there exist a similar end point. Pick different end points")
