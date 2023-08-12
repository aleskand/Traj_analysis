from packages.traj import *
import plotly.graph_objs as go
from plotly_resampler import FigureResampler


class Plot:
    def __init__(self, trajectory, plot_name, plot_scope, debug):
        self.frame_list = trajectory.frame_list
        self.knot_dict = trajectory.knot_dict
        self.untied_list = []
        self.plot_dict = {}
        self.trajectory = trajectory
        self.plot_name = plot_name
        self.plot_scope = plot_scope
        self.debug = debug

    def draw_plot(self):
        """
        Function prepares and processes data in order to draw the plot:
            1. Constructs the plot_dict, where the frames around which the knot was tied are keys and the values are
               the knot core ranges in successive frames. The knot core ranges are necessary to draw the plot.
            2. Draws the plot of knot core range in the whole trajectory.

        Returns:
            The plot.
        """
        self.plot_dict = self.prepare_data_to_plot()
        self.generate_plot()

    def prepare_data_to_plot(self):
        """
        Function creates a plot dict which is necessary to draw the plot. The keys are the frames around which the
        knot was tied, the values are the knot core ranges in successive frames. After knotting, the function
        calculates the knot value for 10 frames every one, then every plot_scope till the frame, where the knot
        was unknotted.

        Returns: plot_dict
        """
        self.plot_dict = {}

        # a variable that determines whether the knot core for the last frame should be counted
        if_end = False
        # a variable that determines whether the knot core for the last frame has been counted
        end_plot = True

        for i, frame in enumerate(self.knot_dict):
            self.plot_dict[frame] = []
            end = self.knot_dict[frame][1]
            if end is None:
                end = self.trajectory.max_frame
                if_end = True

            # calculating the first 10 frames every 1
            for j in range(frame, frame + 11):
                self.plot_dict[frame].append(knotcore_len(j, self.trajectory.lx, self.trajectory.closure,
                                                          self.trajectory.tries, self.trajectory.max_cross))

            # calculating remaining frames every plot_scope
            for j in range(frame + 10 + 100, end, self.plot_scope):
                self.plot_dict[frame].append(knotcore_len(j, self.trajectory.lx, self.trajectory.closure,
                                                          self.trajectory.tries, self.trajectory.max_cross))
                if j == self.trajectory.max_frame:
                    end_plot = False

        if if_end and end_plot:
            self.plot_dict[self.frame_list[-1]].append(knotcore_len(self.trajectory.max_frame, self.trajectory.lx,
                                                                    self.trajectory.closure, self.trajectory.tries,
                                                                    self.trajectory.max_cross))
        return self.plot_dict

    def generate_plot(self):
        """
        Function generates the plot in the following steps. It constructs lists of values for x and two lines y based
        on plot_dict. Then, it gathers information about the knot's tying method and loop position and subsequently
        plots (by using helper functions) the segments on the graph.
        """

        def draw_plot_section(name, x_range, y_range, mode, line, marker, showlegend, hover, fill, fillcolor):
            """
            Auxiliary function, which draws a fragment of the plot.
            name (str):
                    Name of the plotting line. I will be displayed on the plot.
            x_range (list, [int, int]):
                    X-axis data range.
            y_range (list, [int, int]):
                    Y-axis data range.
            mode (str):
                    From plotly.graph_objs._scatter.Scatter documentation:
                    Determines the drawing mode for this scatter trace. If the provided `mode` includes "text" then
                    the `text` elements appear at the coordinates. Otherwise, the `text` elements appear on hover.
                    If there are less than 20 points and the trace is not stacked then the default is "lines+markers".
                    Otherwise, "lines".
            line (dict):
                    From plotly.graph_objs._scatter.Scatter documentation:
                    :class:`plotly.graph_objects.scatter.Line` instance or dict with compatible properties.
            marker (dict):
                    From plotly.graph_objs._scatter.Scatter documentation:
                    :class:`plotly.graph_objects.scatter.Marker` instance or dict with compatible properties.
            showlegend (bool):
                    Determines whether an item corresponding to this trace is shown in the legend.
            hover (str):
                    Additional info, which will be displayed after hovering over the point of the plotting line
                    with the cursor.
            fill (str):
                    Sets the area to fill with a solid color. For specific info about possible parameters look at
                    plotly.graph_objs._scatter.Scatter documentation (args 'fill').
            fillcolor (str):
                     Sets the fill color.
            """
            return fig.add_trace(go.Scatter(
                name=name,
                x=x_range,
                y=y_range,
                mode=mode,
                line=line,
                marker=marker,
                showlegend=showlegend,
                hovertemplate=hover,
                fill=fill,
                fillcolor=fillcolor))

        def draw_legend_element(mode, name, l_color, dash):
            """
            Auxiliary function, which draws a fragment of chart legend.
            mode (str):
                    From plotly.graph_objs._scatter.Scatter documentation:
                    Determines the drawing mode for this scatter trace. If the provided `mode` includes "text" then
                    the `text` elements appear at the coordinates. Otherwise, the `text` elements appear on hover.
                    If there are less than 20 points and the trace is not stacked then the default is "lines+markers".
                    Otherwise, "lines".
            name (str):
                    Name of the element displayed in the legend.
            l_color (str):
                    Color of the displayed line.
            dash (str):
                    Type of the chart line.
            """
            return fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode=mode,
                name=name,
                line=dict(color=l_color, width=1, dash=dash)))

        # dictionary of the colors use in plot
        color_dict = {'3_1': ['rgba(0,128,0, 0.6)', 'green'],
                      '4_1': ['rgba(255,165,0, 0.8)', 'orange'],
                      '5_2': ['rgba(148,0,211, 0.6)', 'darkviolet'],
                      '6_1': ['rgba(100,149,237, 0.6)', 'cornflowerblue'],
                      }

        # dictionary of the colors use in first 10 frames of every knotted plot
        color_dict10 = {'3_1': ['rgba(0,128,0, 0.4)', 'green'],
                        '4_1': ['rgba(255,165,0, 0.6)', 'orange'],
                        '5_2': ['rgba(148,0,211, 0.4)', 'darkviolet'],
                        '6_1': ['rgba(100,149,237, 0.4)', 'cornflowerblue'],
                        }

        fig = (go.Figure())
        legend_entries = {}
        x = []
        i = 0

        # without it, the plot "draws badly"
        fig.add_trace(go.Scatter(
            name=".",
            x=[x for x in range(0, self.trajectory.max_frame, 1)],
            y=[0] * self.trajectory.max_frame,
            mode='lines',
            line=dict(color='white', width=0.5),
            hoverinfo='none',
            showlegend=False))

        # information of the protein length
        fig.update_layout(
            annotations=[
                go.layout.Annotation(
                    text="Protein length: " + str(self.trajectory.prot_len),
                    xref="paper",
                    yref="y",
                    x=1.09,
                    y=self.trajectory.prot_len,
                    showarrow=False,
                    font=dict(color="black", size=15))])

        while i <= self.trajectory.max_frame:
            x.append(i)

            if i in self.plot_dict:
                fig.add_trace(go.Scatter(
                    name='',
                    x=x,
                    y=[None] * len(x),
                    marker=dict(color="#444"),
                    showlegend=False))

                y_upper = []
                y_lower = []
                temp_x = []
                j = 0

                for tup in self.plot_dict[i]:
                    if tup is not None and tup != 0:
                        y_lower.append(tup[0])
                        y_upper.append(tup[1])
                    else:
                        y_lower.append(None)
                        y_upper.append(None)
                    temp_x.append(i + j)

                    if j < 10:
                        j += 1

                    else:
                        if i + j + self.plot_scope < self.trajectory.max_frame:
                            j += self.plot_scope
                        else:
                            j = self.trajectory.max_frame - i

                knot = self.knot_dict[i][0]

                if knot not in legend_entries:
                    legend_entries[knot] = color_dict[knot][1]

                color = color_dict[knot][0]
                color10 = color_dict10[knot][0]

                if self.knot_dict[i][2] == 1:
                    if self.trajectory.nterminus:
                        # Slipknot
                        draw_plot_section("Slipknot", [i - 1, i], [0, self.plot_dict[i][0][0]],
                                          'lines', dict(dash='dot', color='black', width=1),
                                          dict(color='black', size=6), False, "N-terminus", None, None)
                    else:
                        draw_plot_section("Slipknot", [i - 1, i], [0, self.plot_dict[i][0][1]],
                                          'lines', dict(dash='dot', color='black', width=1),
                                          dict(color='black', size=6), False, "N-terminus", None, None)

                else:
                    # Normally
                    if self.trajectory.nterminus:
                        draw_plot_section("Normally", [i - 1, i], [0, self.plot_dict[i][0][0]],
                                          'lines', dict(dash='dash', color='black', width=1),
                                          dict(color='black', size=6), False, "N-terminus", None, None)
                    else:
                        draw_plot_section("Normally", [i - 1, i], [0, self.plot_dict[i][0][1]],
                                          'lines', dict(dash='dash', color='black', width=1),
                                          dict(color='black', size=6), False, "C-terminus", None, None)

                fig.update_layout(hoverlabel=dict(font_size=14))

                info = ''
                loop_color_l = "black"
                loop_color_u = "black"

                if len(self.knot_dict[i]) == 5:
                    # behavior of the loop
                    if self.knot_dict[i][4] == 0:
                        if self.trajectory.nterminus:
                            loop_color_u = "blue"
                        else:
                            loop_color_l = "blue"
                        info = "loop tightens"
                    if self.knot_dict[i][4] == 1:
                        if self.trajectory.nterminus:
                            loop_color_u = "green"
                        else:
                            loop_color_l = "green"
                        info = "loop in place"
                    if self.knot_dict[i][4] == 2:
                        if self.trajectory.nterminus:
                            loop_color_u = "red"
                        else:
                            loop_color_l = "red"
                        info = "loop expands"

                # plot rest of the plot every scope_plot
                draw_plot_section(knot, temp_x[10:], y_lower[10:], 'lines+markers', dict(color='black', width=1),
                                  dict(symbol='circle', size=6), False, None, 'tozeroy', 'rgba(0,0,0,0)')

                draw_plot_section(knot, temp_x[10:], y_upper[10:], 'lines+markers', dict(color='black', width=1),
                                  dict(symbol='circle', size=6), False, None, 'tonexty', color)

                # draw first 10 frames
                draw_plot_section(knot, temp_x[:11], y_lower[:11], 'lines+markers', dict(color=loop_color_l, width=1),
                                  dict(symbol='circle', size=6), False, None, 'tozeroy', 'rgba(0,0,0,0)')

                draw_plot_section(knot, temp_x[:11], y_upper[:11], 'lines+markers', dict(color=loop_color_u, width=1),
                                  dict(symbol='circle', size=6), False, "%{y}, " + info, 'tonexty', color10)

                x = [temp_x[-1] + 1]

            if i == self.trajectory.max_frame:
                fig.add_trace(go.Scatter(
                    name='',
                    x=x,
                    y=[None] * len(x),
                    marker=dict(color="#444"),
                    showlegend=False))
            i += 1

        fig.update_layout(
            yaxis=dict(
                range=[0, self.trajectory.prot_len],
                title='Residue index',
                title_font=dict(
                    color='black',
                    size=15,
                    family='Arial'),
                tickmode='linear',
                dtick=20,
                showticklabels=True,
                ticks='outside',
                ticklen=10,
                tickwidth=2,
                tickcolor='black'),

            xaxis=dict(title='Frame',
                       title_font=dict(color='black', size=15, family='Arial')),
            title='The knot core plot',
            title_font=dict(
                color='black',
                size=20,
                family='Arial'),
            hovermode="x",
            autosize=True)

        fig.update_layout(
            xaxis=dict(
                rangeselector=dict(
                    buttons=list([dict(count=1, stepmode="backward")])),
                rangeslider=dict(visible=True),
                tickformat=',d'))

        fig.update_layout(
            legend=dict(
                y=0.9,
                bordercolor="black",
                borderwidth=1,
                bgcolor="white",
                font=dict(color='black', size=14),
                itemsizing="constant"))

        # Plotting legend
        draw_legend_element(None, 'Loop tightens', 'blue', 'solid')
        draw_legend_element(None, 'Loop in place', 'green', 'solid')
        draw_legend_element(None, 'Loop expands', 'red', 'solid')
        draw_legend_element(None, 'Loop behavior:', 'white', None)
        draw_legend_element(None, "", 'white', None)
        draw_legend_element('lines', 'Normally', 'black', 'dash')
        draw_legend_element('lines', 'Via slipknot', 'black', 'dot')
        draw_legend_element(None, 'Way of looping: ', 'white', None)
        draw_legend_element(None, "", 'white', None)

        for knot in legend_entries:
            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                name=knot,
                marker=dict(size=15, color=legend_entries[knot], symbol='square')))

        draw_legend_element(None, 'Knot type:', 'white', None)

        fig.update_layout(legend=dict(itemsizing='trace'))
        fig.update_traces(line=dict(width=3))

        resampled_fig = FigureResampler(fig)

        resampled_fig.show()
        resampled_fig.write_html(self.plot_name + '.html')
