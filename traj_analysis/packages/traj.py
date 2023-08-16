import math
from packages.knotcore import *
import os


def search_for_the_type_change(start, end, iteration, lx, closure, max_cross, tries, loop):
    """
    Function iterates every specified step, searching for the moments when the type of knot in the trajectory changes.

    Args:
        start (int):
                Number of frame where the search starts.
        end (int):
                Number of frame where the search ends.
        iteration (int: 1, 10 or 100):
                Step every which we perform a trajectory search.
        loop (bool):
                True if looking for looping moment.

    Returns: list of frames, in which a knot type change was detected
    """
    frame_list = []
    if loop:
        knot = "0_1"
    else:
        knot = " "

    for i in range(start, end, iteration):
        kn = knot_type(i, lx, closure, max_cross, tries)

        if loop and str(kn) != knot and str(kn) != '0_1':
            frame_list.append(i)
            if iteration == 10 or iteration == 1:
                break
        if not loop and str(kn) != knot and str(kn) == '0_1' and i != 0:
            frame_list.append(i)
            if iteration == 10 or iteration == 1:
                break
        elif iteration == 10 and i == end - 10:
            frame_list.append(end)
        knot = str(kn)

    return frame_list


def knot_type(i, lx, closure, max_cross, tries):
    """
    Function calculates the Alexander polynomial of the given structure.

    Returns: topology type.
    """
    if closure == 1:
        return alexander([[x, y, z] for x, y, z in lx[i][::2]], closure=closure, run_parallel=False,
                         max_cross=max_cross)
    else:
        res = alexander([[x, y, z] for x, y, z in lx[i][::2]], closure=closure, tries=tries, run_parallel=False,
                        max_cross=max_cross)
        max_value = max(res.values())
        max_keys = [key for key, value in res.items() if value == max_value]
        if len(max_keys) > 1:
            return '0_1'
        else:
            return max_keys[0]


def knotcore_len(i, lx, closure, tries, max_cross):
    """
    Function creates a nxyz file and calculates knot core value in the given frame.

    Returns: knot core value
             None, if the knot core function returns invalid value.
    """

    j = 0
    # creating nxyz file
    with open("frames_temp.nxyz", "a") as output:
        for line in lx[i]:
            output.write(str(j))
            output.write(' ')
            for x in line:
                output.write(' ')
                output.write(str(x))
            output.write('\n')
            j += 1
    output.close()

    # calculate knot core
    knotcore_res = count_knotcore('frames_temp.nxyz', closure=closure, tries=tries, max_cross=max_cross)
    if knotcore_res is None:
        knotcore_res = 0
    else:
        if knotcore_res[1] - knotcore_res[0] == 0:
            knotcore_res = None

    os.remove("frames_temp.nxyz")

    return knotcore_res


def check_after_knotting(start, end, lx, closure, max_cross, tries):
    """
    Function checks if knot is tied on the correct number of frames.

    Returns: 1 if knotted correctly,
            frame, in which the unknot is found otherwise.
    """
    for i in range(start, end):
        kn = knot_type(i, lx, closure, max_cross, tries)
        if kn == '0_1':
            return i
    return 1


def check_knotting(start, end, pc, lx, min_gap, closure, max_cross, tries):
    """
    Function checks if knot is not tied on the given percentage (pc) of min_gap frames. It can be used to check, if
    there were enough frames without a knot before the moment of knotting. Or to check, if the knot was really
    unknotted, or the unknot frames appeared by chance after the consider moment of unknotting.

    Returns: True if not knotted, False otherwise.
    """
    case = math.floor((1-pc) * min_gap)

    for i in range(start, end):
        kn = knot_type(i, lx, closure, max_cross, tries)
        if kn != '0_1' and case < 0:
            return False
        if kn != '0_1':
            case -= 1
    return True


class Traj:
    def __init__(self, lx, prot_len, max_frame, min_gap, scope, min_knot, nterminus, nat_knotcore, closure, tries,
                 max_cross, debug):
        self.lx = lx
        self.prot_len = prot_len
        # maximum tail length for slipknot classification, 2 thresholds for small (below 100 nucleotides) and
        # large (greater than 100 nucleotides) structures
        if self.prot_len < 100:
            self.SLIPKNOT_SIZE = 5
        else:
            self.SLIPKNOT_SIZE = 10
        self.max_frame = max_frame
        self.min_gap = min_gap
        self.scope = scope
        self.min_knot = min_knot
        self.nterminus = nterminus
        self.nat_knotcore = nat_knotcore
        self.closure = closure
        self.tries = tries
        self.max_cross = max_cross
        self.debug = debug
        self.frame_list = []
        self.knot_dict = {}
        self.untied_list = []

    def calculate(self, full_output):
        """
        The function performs the entire analysis process by calling subsequent functions one by one.
        Analysis steps:
            1. Creating frame_list, which holds frames where is possible that the knot was tied.
            2. Creating knot_dict - dictionary, where the results of analysis will be stored. During this process,
               checks are made to find the exact moment of knotting.
            3. Calculating the knot core range for the frames, in which the knot was tied.
            4. Evaluating the way of looping and the behavior of the loop, after knotting.
        All the above information are inserted to the knot_dict and returned.

        Returns:
            Dictionary with the results of the analysis.
        """
        self.frame_list = self.searched_structure(True)
        self.untied_list = self.searched_structure(False)
        self.check_untied_list()
        if self.debug:
            print("Result of first iteration of searching for the possible moments of knotting: ", self.frame_list)
        if len(self.frame_list) != 0:
            self.knot_dict = self.construct_knotdict()
            self.check_knot()
            self.calculate_knotcore()
            self.specify_knotting_style()

            if full_output:
                full_knot_dict = {}
                for key in self.knot_dict:
                    full_knot_dict[key] = {}
                    full_knot_dict[key]["Knot type"] = self.knot_dict[key][0]
                    full_knot_dict[key]["Unknotting frame"] = self.knot_dict[key][1]
                    full_knot_dict[key]["Knot core range"] = self.knot_dict[key][2]
                    if self.knot_dict[key][3] == 0:
                        full_knot_dict[key]["Knotting via slipknot"] = True
                    else:
                        full_knot_dict[key]["Knotting via slipknot"] = False
                    if self.nat_knotcore is not None:
                        if self.knot_dict[key][4] == 0:
                            full_knot_dict[key]["Loop behavior"] = "loop tightens"
                        if self.knot_dict[key][4] == 1:
                            full_knot_dict[key]["Loop behavior"] = "loop is in place"
                        if self.knot_dict[key][4] == 2:
                            full_knot_dict[key]["Loop behavior"] = "loop expands"
                    else:
                        full_knot_dict[key]["Loop behavior"] = "No rating. The knot core range the native form of the" \
                                                               " structure was not given."
                return full_knot_dict
            else:
                return self.knot_dict
        else:
            return None

    def searched_structure(self, knotting):
        """
        Function searches the trajectory to find the moment of change from unknot to knot or from knot to unknot.
        First, it searches every 100 frames, then every 10 over the 100 frames before the frame, which was found in
        the previous iteration, then every 1 frame.
        Thanks to this, the function finds possible moments of knotting or unknotting, which will be carefully analyzed
        further on. The function by default ignores the possible momentary creation of the knot (for less than 100
        frames), because its purpose is to find those moments when a stable knot is created.

        Args:
            knotting (bool):
                    True, if looking dor the moments of knotting.
                    False, if looking for thr moments of unknotting.

        Returns: list of frames, where the knot is likely to have tied or untied.
        """
        # searching every 100 frames
        frame_list_100 = search_for_the_type_change(0, len(self.lx), 100, self.lx, self.closure, self.max_cross,
                                                    self.tries, knotting)

        # searching every 10 frames
        frame_list_10 = []
        for j in range(len(frame_list_100)):
            frame = search_for_the_type_change(frame_list_100[j] - 90, frame_list_100[j]-10, 10, self.lx,
                                               self.closure, self.max_cross, self.tries, knotting)
            if len(frame) == 0:
                frame_list_10.append(frame_list_100[j])
            else:
                frame_list_10.append(frame[0])

        # searching every 1 frame
        frame_list_1 = []
        for j in range(len(frame_list_10)):
            frame = search_for_the_type_change(frame_list_10[j] - 9, frame_list_10[j]-1, 1, self.lx,
                                               self.closure, self.max_cross, self.tries, knotting)
            if len(frame) == 0:
                frame_list_1.append(frame_list_10[j])
            else:
                frame_list_1.append(frame[0])

        return frame_list_1

    def check_untied_list(self):
        """
        Function checks the untied list. If the length of the 'untied' list is not equal to the length of the
        'frame_list', it means that either the knot was tied from the beginning, in which case a 0 is appended to the
        'frame_list', or the knot remained unresolved until the end of the trajectory. In such a situation, the value
        None is inserted into the 'untied' list. This information will be used by subsequent functions to indicate that
        the last tied knot, after being tied, remained present until the end of the trajectory.
        """
        if len(self.untied_list) != len(self.frame_list):
            if len(self.frame_list) == 1 and len(self.untied_list) == 0:
                self.untied_list.append(None)
            else:
                if self.untied_list[0] < self.frame_list[0]:
                    self.frame_list.insert(0, 0)
                if self.frame_list[-1] > self.untied_list[-1]:
                    self.untied_list.append(None)

    def construct_knotdict(self):
        """
        Function analyzes the data and calculates the knot_dict, by performing the following operations:
        function looks for a stable knot formation in the vicinity of frames from the "frame_list", by checking whether
        there are at least "min_gap" frames that are 80% unknotted before the consider frame of knotting. If this
        condition is met, it checks whether there are at least "scope" frames that are knotted after the consider
        knotting frame. If the selected frame passes this verification successfully, function considers that a stable
        knot has formed in it, and packages it to the dictionary.

        If the second condition is not met, the function starts the next iteration of checking from the next frame
        after the one in which an unknot was detected. The verification process is conducted following the same rules
        as before, with the exception that it no longer checks the condition for untangling before the knotting moment.
        There can be maximum 10 such iterations. This limitation is introduced to prevent the function from getting
        stuck in a loop.

        If it turns out that in no frame a stable knot was formed, lasting for at least 10 frames, then function
        concludes, that there is no knot in the analyzed vicinity, which was initially considered a possible knotting
        location. The index of this frame is added to auxiliary_list. This list holds the indexes, of frames in which
        the analyzed knot is unstable. Later, the corresponding to untie moments for those frames will be removed from
        the 'untied_list', thus eliminating all non-stable looping moments from the analysis of knot formation.

        The last step of this function is to check whether the frames in which the knot untied itself meet the
        specified condition (whether at least 50% of the frames within 'CHECK_LEN' (default=10) are untied after the
        moment of unknotting).

        You can change the percentage value:
                PC_KNOTTING = 0.8 - The percentage of frames within 'min_gap' that need to be untied before knot
                                    formation to consider that the knot has actually formed.
                PC_UNKNOTTING = 0.5 - The percentage of frames within 'CHECK_LEN' that must be untied after the
                                      analyzed untie moment, in order to actually consider it as an untie moment.
                CHECK_LEN = 10 - The number of frames that will be taken into account when checking whether the knot is
                                 resolved.

        Returns: knot_dict - Dictionary of frames in which knot is knotted as keys and a topology type in a lit as
        values. Later on the list will be updated with the result of further analysis.
        """
        knot_dict = {}
        auxiliary_list = []
        PC_KNOTTING = 0.8
        PC_UNKNOTTING = 0.5
        CHECK_LEN = 10
        for frame_index, frame in enumerate(self.frame_list):
            found = False
            # check if there was no knot before the found frame
            if check_knotting(frame - self.min_gap, frame - 1, PC_KNOTTING, self.lx, self.min_gap, self.closure,
                              self.max_cross, self.tries):
                # check if knot is tied on the correct number of frames
                result = check_after_knotting(frame + 1, frame + self.scope, self.lx, self.closure, self.max_cross,
                                              self.tries)
                if result == 1:
                    knot_dict[frame] = []
                else:
                    # knot is not tied correctly, further checks, but maximum 10 times
                    counter = 0
                    while counter < 10:
                        kn = knot_type(result + 1, self.lx, self.closure, self.max_cross, self.tries)
                        if kn != '0_1':
                            check = check_after_knotting(result + 2, result + self.scope - 1, self.lx, self.closure,
                                                         self.max_cross, self.tries)
                            if check == 1:
                                # knot find in this frame is correct
                                knot_dict[result + 1] = []
                                found = True
                                break
                            else:
                                # knot is not tied correctly, further checks,
                                result = check
                                counter += 1
                        else:
                            # knot in the frame nr result+1 is an unknot
                            counter += 1
                            result += 2
                    if not found:
                        auxiliary_list.append(frame_index)
            else:
                auxiliary_list.append(frame_index)

        # calculating knot type
        for frame in knot_dict:
            kn = knot_type(frame, self.lx, self.closure, self.max_cross, self.tries)
            knot_dict[frame] = [kn]

        # updating the untied_list
        for f in auxiliary_list:
            if 0 <= f < len(self.untied_list):
                del self.untied_list[f]

        # checking and inserting to knot_dict the frame number where the knot untied
        for i, un_frame in enumerate(self.untied_list):
            if un_frame is None:
                break
            if check_knotting(un_frame + 1, un_frame + CHECK_LEN, PC_UNKNOTTING, self.lx, CHECK_LEN, self.closure,
                              self.max_cross, self.tries):
                frame = list(knot_dict.keys())[i]
                knot_dict[frame].append(un_frame)
            else:
                find = True
                next_frame = un_frame + 11
                while find and next_frame + 10 < self.max_frame:
                    if check_knotting(next_frame, next_frame + CHECK_LEN, PC_UNKNOTTING, self.lx, CHECK_LEN,
                                      self.closure, self.max_cross, self.tries):
                        frame = list(knot_dict.keys())[i]
                        knot_dict[frame].append(un_frame)
                        find = False
                    else:
                        next_frame += 11

        # Checking whether the last knot has a recorded untie moment or if it is tied until the end of the
        # and needs this information to be added.
        last_frame = list(knot_dict.keys())[-1]
        if len(knot_dict[last_frame]) == 1:
            knot_dict[last_frame].append(None)

        return knot_dict

    def check_knot(self):
        """
        Function checks if the distance between the formation of the knot and its unknotting frame meets the condition
        adopted in the analysis, i.e. whether it is greater than min_knot. If not, it modifies the knot_dict.
        """
        too_short = []
        for frame in self.knot_dict.keys():
            end = self.knot_dict[frame][1]
            if end is None:
                end = self.max_frame
            else:
                end = int(self.knot_dict[frame][1])
            if end - int(frame) < self.min_knot:
                too_short.append(frame)

        for nr in too_short:
            del self.knot_dict[nr]

    def calculate_knotcore(self):
        """
        The function calculate the knot core range in frames, where the knot was tied and inserts the results into
        the knot_dict, at the same time validating the knot core range. If it receives an incorrect value, it searches
        for the correct value in the next 10 frames. Incorrect values are written to the list 'keys_to_modify' and
        corrected at the end of the function.
        """
        keys_to_modify = []
        for frame in self.knot_dict:
            er = False
            knotcore = knotcore_len(frame, self.lx, self.closure, self.tries, self.max_cross)
            try:
                if isinstance(knotcore, int):
                    raise TypeError("Knot core value can not be 0.")
                elif knotcore is None:
                    raise ValueError("Knot core length must be greater than 6.")
            except (TypeError, ValueError) as e:
                if self.debug:
                    print(f"Exception occurred: {e}. Analysis of frame {frame} is processed. Looking for the knot in "
                          f"successive frames.")
                er = True
            if er:
                # Invalid knot core in frame, looking for the correct value in subsequent frames, but maximum in 10
                # frames
                for i in range(frame + 1, frame + 10):
                    knotcore = knotcore_len(i, self.lx, self.closure, self.tries, self.max_cross)
                    if type(knotcore) is tuple:
                        if knotcore[1] - knotcore[0] > 6:
                            keys_to_modify.append((frame, i, knotcore))
                            break
            else:
                self.knot_dict[frame].append(knotcore)

        # Modify the dictionary
        for old_frame, new_frame, knotcore_value in keys_to_modify:
            value = self.knot_dict[old_frame]
            self.knot_dict[new_frame] = value
            self.knot_dict[new_frame].append(knotcore_value)
            del self.knot_dict[old_frame]

        self.knot_dict = dict(sorted(self.knot_dict.items()))

    def specify_knotting_style(self):
        """
        The function evaluates the way of knotting, the behavior of the loop and inserts the results into the knot_dict.

        Returns:
            1. Way of knotting:
                0, if the loop forms via slipknot.
                1, if the loop forms directly.
            2. Behavior of the loop:
                0, if the loop tightens.
                1, if the loop in place.
                2, if the loop expands.
        """

        for frame in self.knot_dict:
            knotcore = self.knot_dict[frame][2]

            # rating the way of knotting
            if self.nterminus:
                if knotcore[0] >= self.SLIPKNOT_SIZE:
                    # N-terminus, slipknot
                    self.knot_dict[frame].append(0)
                else:
                    # N-terminus, directly
                    self.knot_dict[frame].append(1)
            else:
                if self.prot_len - knotcore[1] >= self.SLIPKNOT_SIZE:
                    # C-terminus, slipknot
                    self.knot_dict[frame].append(0)
                else:
                    # C-terminus, directly
                    self.knot_dict[frame].append(1)

            if self.nat_knotcore is not None:
                if self.nterminus:
                    diff = (self.prot_len - knotcore[1]) / (self.prot_len - self.nat_knotcore[1])
                else:
                    diff = knotcore[0] / self.nat_knotcore[0]

                # rating the behavior of the loop
                if diff < 0.5:
                    # loop is tightens
                    self.knot_dict[frame].append(0)
                elif 0.5 <= diff <= 1.5:
                    # loop in place
                    self.knot_dict[frame].append(1)
                elif diff > 1.5:
                    # loop expands
                    self.knot_dict[frame].append(2)

        if self.debug and self.nat_knotcore is None:
            print("The knot core range of the native form of the structure was not given. Function didn't rate"
                  " the behavior of the loop.")
