import h5py
import numpy as np
import pandas as pd
import multiprocessing as mp
from functools import partial
import traceback
import scipy.signal as sg
from scipy.stats import zscore
from scipy.stats import _stats_py
from MEA_analysis import RF_code_new as rf
from scipy.spatial import KDTree
from scipy import signal, ndimage
import networkx as nx
import plotly.express as px
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
import plotly.graph_objects as go
import re
import random
from ipywidgets import Button, HBox, VBox


# def zmap(scores, compare, axis=0, ddof=0):
#     """
#     Calculate the relative z-scores.
#
#     Return an array of z-scores, i.e., scores that are standardized to
#     zero mean and unit variance, where mean and variance are calculated
#     from the comparison array.
#
#     Parameters
#     ----------
#     scores : array_like
#         The input for which z-scores are calculated.
#     compare : array_like
#         The input from which the mean and standard deviation of the
#         normalization are taken; assumed to have the same dimension as
#         `scores`.
#     axis : int or None, optional
#         Axis over which mean and variance of `compare` are calculated.
#         Default is 0. If None, compute over the whole array `scores`.
#     ddof : int, optional
#         Degrees of freedom correction in the calculation of the
#         standard deviation. Default is 0.
#
#
#     Returns
#     -------
#     zscore : array_like
#         Z-scores, in the same shape as `scores`.
#
#     Notes
#     -----
#     This function preserves ndarray subclasses, and works also with
#     matrices and masked arrays (it uses `asanyarray` instead of
#     `asarray` for parameters).
#
#     Examples
#     --------
#     >>> from scipy.stats import zmap
#     >>> a = [0.5, 2.0, 2.5, 3]
#     >>> b = [0, 1, 2, 3, 4]
#     >>> zmap(a, b)
#     array([-1.06066017,  0.        ,  0.35355339,  0.70710678])
#
#     """
#     a = np.asanyarray(compare)
#
#     if a.size == 0:
#         return np.empty(a.shape)
#
#     mn = np.median(a, axis=axis, keepdims=True)
#     std = a.std(axis=axis, ddof=ddof, keepdims=True)
#     if axis is None:
#         isconst = (a.item(0) == a).all()
#     else:
#         isconst = (_stats_py._first(a, axis) == a).all(axis=axis, keepdims=True)
#
#     # Set std deviations that are 0 to 1 to avoid division by 0.
#     std[isconst] = 1.0
#     z = (scores - mn) / std
#     # Set the outputs associated with a constant input to nan.
#     z[np.broadcast_to(isconst, z.shape)] = np.nan
#     return z


def create_trigger_signal(file, nr_frames, stimbegin=True):
    trigger = np.arange(0, get_nr_of_frames(file), 1)
    if stimbegin:
        trigger = trigger + find_stimbegin(file, nr_frames)
    return trigger


def find_stimbegin(file, nr_frames):
    trace = brw_import_chunks(file, 0, nr_frames)
    half_voltage = np.max(trace[1, :]) * 0.7
    print(half_voltage)
    channel_log = trace[1, :] > half_voltage
    peaks = sg.find_peaks(channel_log, height=1, plateau_size=2)
    peaks_left = peaks[1]["left_edges"]
    return peaks_left[0]


def get_nr_of_frames(file):
    with h5py.File(file, "r") as f:
        rec_vars = f.require_group("3BRecInfo/3BRecVars/")
        n_frames_total = rec_vars["NRecFrames"][()][0]
    return n_frames_total


def brw_import_chunks(file, f0, f1):
    nr_fr = f1 - f0
    with h5py.File(file, "r") as f:
        rec_vars = f.require_group("3BRecInfo/3BRecVars/")
        n_frames_total = rec_vars["NRecFrames"][()][0]
        n_rec_ch = int(1.0 * f["3BData/Raw"].shape[0] / n_frames_total)
        start_frame = f0 * n_rec_ch
        end_frame = f1 * n_rec_ch

        data = f["/3BData/Raw"][start_frame:end_frame]
        data = np.reshape(data, (n_rec_ch, nr_fr), order="F").astype(int)
        return data


def brw_import_chunks_loop(file, f0, f1):
    nr_fr = f1 - f0
    rec_vars = file.require_group("3BRecInfo/3BRecVars/")
    n_frames_total = rec_vars["NRecFrames"][()][0]
    n_rec_ch = int(1.0 * file["3BData/Raw"].shape[0] / n_frames_total)
    start_frame = f0 * n_rec_ch
    end_frame = f1 * n_rec_ch
    data = np.zeros((int(end_frame - start_frame)), dtype=int)
    data[:] = file["/3BData/Raw"][start_frame:end_frame]
    data = np.reshape(data, (n_rec_ch, nr_fr), order="F")
    return data


def sta_main_loop(spikes, trigger, file, nr_pre_bins=20, nr_post_bins=5, nr_boxes=1):
    h5f = h5py.File(file, "r")
    bins = nr_pre_bins + nr_post_bins
    # Creates fake spikes_counts data all = 1

    sta = np.zeros((nr_boxes, bins), dtype=int)

    current_trigger = nr_pre_bins  # The trigger to look for a spike
    left_trigger = 0  # last trigger before the spike to be considered
    right_trigger = bins  # the last trigger after the spike to be considered

    for spike in spikes:
        # print(spike)

        # Check if the spiketime is smaller than the trigger, go to next spike
        # (this should really just happen for the first few spikes)
        if right_trigger > np.shape(trigger)[0]:
            print("break")
            break

        if spike < trigger[current_trigger]:
            print("smaller")
            continue

        # This checks if the STA bins would run over the stimulus duration in the end
        # and breaks the loop in that case

        # If the spike is spiketime is larger than the trigger, a while loop starts.
        # This while loop basically checks of the spike is smaller than the next trigger
        # (ans so lies between trigger and next trigger) or if its even larger and the
        # a larger current trigger has to be considered

        # if spike > trigger[current_trigger]:
        #     try:
        #         while spike > trigger[current_trigger + 1]:
        #             current_trigger = current_trigger + 1
        #             left_trigger = left_trigger + 1
        #             right_trigger = right_trigger + 1
        #             print("stuck")
        #     except IndexError:
        #         print("index", spike, current_trigger)
        #         break
        left_trigger = spike - nr_pre_bins
        right_trigger = spike + nr_post_bins
        # When the while loop has ended, it is sure that the spike happend between the
        # current trigger and the next trigger and we can allocate the STA bins from
        # the noise sequence

        # In this case the sequence will be multiplied with the number of unique spikes

        # within this bin. This step is absent from the "normal" version of the loop
        if right_trigger - left_trigger < bins - 1:
            print("bins")
            left_trigger = left_trigger - (right_trigger - left_trigger)

        try:
            # "Normal version"
            trace = brw_import_chunks_loop(h5f, left_trigger, right_trigger)
            # trace_norm = zscore(trace, axis=1, nan_policy="omit")
            detrended = signal.detrend(trace, bp=np.arange(0, 50, 50), axis=1)
            win = signal.windows.gaussian(2, std=1)

            filtered = ndimage.convolve1d(
                detrended, win, axis=0, mode="constant"
            ) / sum(win)

            # trace_norm = zscore(trace, axis=1, nan_policy="omit")
            # trace = (trace - np.mean(trace, axis=1)[np.newaxis].transpose()) / np.std(
            #     trace, axis=1
            # )[np.newaxis].transpose()
            # if np.any(np.isnan(trace)):
            #     trace[np.isnan(trace)] = 0
            # trace[trace > 2200] = np.mean(trace)
            if np.max(np.absolute(np.diff(filtered))) > 1000:
                continue
            else:
                sta = sta + filtered

        except ValueError:
            print("Value", left_trigger, right_trigger)
    h5f.close()
    return sta


def sta_main_loop_snippets(
    spikes, trigger, file, nr_pre_bins=20, nr_post_bins=5, nr_boxes=1
):
    h5f = h5py.File(file, "r")
    bins = nr_pre_bins + nr_post_bins
    # Creates fake spikes_counts data all = 1
    nr_spikes = spikes.shape[0]
    sta = np.zeros((nr_spikes, nr_boxes, bins), dtype=int)

    current_trigger = nr_pre_bins  # The trigger to look for a spike
    left_trigger = 0  # last trigger before the spike to be considered
    right_trigger = bins  # the last trigger after the spike to be considered

    for spike, spike_idx in zip(spikes, range(nr_spikes)):
        # print(spike)

        # Check if the spiketime is smaller than the trigger, go to next spike
        # (this should really just happen for the first few spikes)
        if right_trigger > np.shape(trigger)[0]:
            print("break")
            break

        if spike < trigger[current_trigger]:
            print("smaller")
            continue

        # This checks if the STA bins would run over the stimulus duration in the end
        # and breaks the loop in that case

        # If the spike is spiketime is larger than the trigger, a while loop starts.
        # This while loop basically checks of the spike is smaller than the next trigger
        # (ans so lies between trigger and next trigger) or if its even larger and the
        # a larger current trigger has to be considered

        # if spike > trigger[current_trigger]:
        #     try:
        #         while spike > trigger[current_trigger + 1]:
        #             current_trigger = current_trigger + 1
        #             left_trigger = left_trigger + 1
        #             right_trigger = right_trigger + 1
        #             print("stuck")
        #     except IndexError:
        #         print("index", spike, current_trigger)
        #         break
        left_trigger = spike - nr_pre_bins
        right_trigger = spike + nr_post_bins
        # When the while loop has ended, it is sure that the spike happend between the
        # current trigger and the next trigger and we can allocate the STA bins from
        # the noise sequence

        # In this case the sequence will be multiplied with the number of unique spikes

        # within this bin. This step is absent from the "normal" version of the loop
        if right_trigger - left_trigger < bins - 1:
            print("bins")
            left_trigger = left_trigger - (right_trigger - left_trigger)

        try:
            # "Normal version"
            trace = brw_import_chunks_loop(h5f, left_trigger, right_trigger)
            # trace_norm = zscore(trace, axis=1, nan_policy="omit")
            detrended = signal.detrend(trace, bp=np.arange(0, 50, 50), axis=1)
            win = signal.windows.gaussian(2, std=1)

            filtered = ndimage.convolve1d(
                detrended, win, axis=0, mode="constant"
            ) / sum(win)

            # trace_norm = zscore(trace, axis=1, nan_policy="omit")
            # trace = (trace - np.mean(trace, axis=1)[np.newaxis].transpose()) / np.std(
            #     trace, axis=1
            # )[np.newaxis].transpose()
            # if np.any(np.isnan(trace)):
            #     trace[np.isnan(trace)] = 0
            # trace[trace > 2200] = np.mean(trace)
            sta[spike_idx, :, :] = filtered

        except ValueError:
            print("Value", left_trigger, right_trigger)
    h5f.close()
    return sta


def find_neuron(
    sta, boxes=64, tries=20, spike_window=10, nr_neighbours=2, return_counts=False
):
    sta[0, :] = np.mean(sta)
    # sta = zscore(sta, axis=1)
    # sta[np.isnan(sta)] = np.nanmean(sta)
    sta_test = sta.copy()

    corr_count = np.ones((boxes * boxes))
    sta_shape = np.shape(sta)
    corr = np.zeros((boxes * boxes, sta_shape[1]))
    corr_neg = np.zeros((boxes * boxes, sta_shape[1]))
    mins = np.zeros(tries)
    tempreference_store = []
    min_loc_store = []
    for i in range(tries):

        min_loc = np.unravel_index(sta_test.argmin(), sta_test.shape)
        min_loc_store.append(min_loc)
        mins[i] = np.min(sta_test)
        # print(min_loc)
        try:
            start = min_loc[1] - spike_window
            if start < 0:
                start = 0
            end = min_loc[1] + spike_window
            if end > sta_shape[1]:
                end = sta_shape[1]
            tempreference = sta[min_loc[0], start:end]
            tempreference_store.append(tempreference)
            neighbours_idx = rf.ident_surrounding_channels(
                min_loc[0], nr_neighbours=nr_neighbours, boxes=boxes
            )
            # neighbours_idx = neighbours_idx[nr_neighbours - 1 : nr_neighbours + 1, :]
            # neighbours_idx[neighbours_idx==min_loc[0]] = np.NaN
            neighbours_idx = neighbours_idx[~np.isnan(neighbours_idx)].astype(int)
            # print(neighbours_idx)
            neighbour_channels = sta[neighbours_idx, :]
            nr_actual_neighbours = np.shape(neighbour_channels)[0]
            corr_temp = np.zeros((nr_actual_neighbours, sta_shape[1], 2))
            for neighbour in range(nr_actual_neighbours):
                corr_temp[neighbour, :, 0] = sg.correlate(
                    neighbour_channels[neighbour, :],
                    tempreference,
                    method="direct",
                    mode="same",
                )
                # print(corr_temp)
            corr_temp[:, :, 1] = corr[neighbours_idx, :]
            corr[neighbours_idx, :] = np.max(corr_temp, axis=2)
            corr_temp[:, :, 1] = corr_neg[neighbours_idx, :]
            corr_neg[neighbours_idx, :] = np.min(corr_temp, axis=2)
            corr_count[neighbours_idx] = corr_count[neighbours_idx] + 1
            sta_test[min_loc[0], :] = 0
            sta_test[:, min_loc[1]] = 0
        except ValueError:
            pass

    if return_counts:
        return corr, corr_neg, corr_count, mins
    else:
        # corr = corr / corr_count[np.newaxis].transpose()
        # corr_neg = corr_neg / corr_count[np.newaxis].transpose()

        corr_max = np.max(corr, axis=1)
        corr_max = corr_max / np.max(corr_max)
        corr_max = np.reshape(corr_max, (boxes, boxes))

        corr_min = np.min(corr_neg, axis=1)
        corr_min = corr_min / np.min(corr_min)
        corr_min = np.reshape(corr_min, (boxes, boxes))
        return corr_max, corr_min, mins, tempreference_store, min_loc_store


def filter_neurons(
    stas,
    nr_filter_steps=2,
    neighbours=2,
    filter_threshold=0.2,
    max_channel_nr=400,
    threshold=0.05,
    filter=False,
    filter_array=None,
):
    nr_neurons = stas.shape[0]
    combined_all = np.zeros((nr_neurons, stas.shape[1]))
    positions = np.empty(nr_neurons, dtype=object)
    channel_x = np.sqrt(stas.shape[1]).astype(int)
    for neuron in range(nr_neurons):
        sta = stas[neuron, :, :].copy()
        positions[neuron], combined_all[neuron, :] = filter_neuron(
            sta,
            filter=filter,
            nr_filter_steps=nr_filter_steps,
            neighbours=neighbours,
            filter_threshold=filter_threshold,
            max_channel_nr=max_channel_nr,
            threshold=threshold,
            filter_array=filter_array,
        )

    combined_all_re = np.reshape(combined_all, (nr_neurons, channel_x, channel_x))
    return positions, combined_all_re


def filter_neuron(sta, filter=False, **kwargs):

    channel_x = np.sqrt(sta.shape[0]).astype(int)
    sta[0:2, :] = 0
    sta_re = np.reshape(sta, (channel_x, channel_x, sta.shape[1]))
    diffs = np.max(np.absolute(np.diff(sta, 1, axis=1)), axis=1)
    diffs = diffs / np.max(diffs)

    maxes = np.max(sta, axis=1)
    maxes = maxes / np.max(maxes)
    mins = np.abs(np.min(sta, axis=1))
    mins = mins / np.max(mins)
    combined = diffs * maxes * mins
    if kwargs["threshold"] == "auto":
        threshold = np.quantile(combined, 0.95)
    else:
        threshold = kwargs["threshold"]
    # combined = zmap(combined, combined)

    combined[combined < threshold] = np.NaN
    # combined[combined<np.quantile(combined, 0.95)] = np.NaN
    for i in range(kwargs["nr_filter_steps"]):
        indices = np.where(~np.isnan(combined))[0]
        indices_penalty = np.empty_like(indices).astype(float)

        for index, i in zip(indices, range(indices.shape[0])):
            s_channels = rf.ident_surrounding_channels(
                index, nr_neighbours=kwargs["neighbours"], boxes=channel_x
            )
            if filter:

                s_channels = s_channels[kwargs["filter_array"]]

            s_channels = s_channels[~np.isnan(s_channels)].astype(int)

            indices_penalty[i] = np.sum(np.isnan(combined[s_channels])).astype(float)
            indices_penalty[i] = 1 - indices_penalty[i] / np.shape(s_channels)[0]

        indices_kick = indices[indices_penalty < kwargs["filter_threshold"]]
        combined[indices_kick] = np.NaN

    if np.sum(~np.isnan(combined)) > kwargs["max_channel_nr"]:
        combined[:] = np.NaN
    reshaped_points = np.reshape(combined, (channel_x, channel_x))
    y_dots, x_dots = np.where(~np.isnan(reshaped_points))
    dots_combined = np.zeros((np.shape(y_dots)[0], 2))
    dots_combined[:, 1] = x_dots
    dots_combined[:, 0] = y_dots
    position = dots_combined.astype(int)

    return position, combined


def trace_neuron(sta_re, position, bins=50, point_radius=2, soma_radius=10):
    sta_re_test = sta_re.copy()
    if bins == "auto":
        bins = position.shape[0]
    # bins = np.shape(sta_re)[2]
    min_loc = np.zeros((bins, 4), dtype=int)
    min_loc_false = np.ones(bins, dtype=bool)
    for i in range(bins):
        if np.min(sta_re_test) < 0:
            min_loc[i, :3] = np.unravel_index(sta_re_test.argmin(), sta_re_test.shape)
            min_loc[i, 3] = np.min(sta_re_test)
            sta_re_test[min_loc[i, 0], min_loc[i, 1], :] = 0
        else:
            min_loc_false[i] = False

    min_loc = min_loc[min_loc_false, :]
    min_loc[:, 2] = min_loc[:, 2] - min_loc[0, 2]
    # kd_tree = KDTree(min_loc[:, :2])
    # pairs = kd_tree.query_pairs(r=point_radius)
    #
    # min_loc_filtered = min_loc[np.unique(np.asarray(list(pairs))), :]
    min_loc_sorted = min_loc[np.argsort(min_loc[:, 2]), :]

    # min_loc_sorted[:, 2] = min_loc_sorted[:, 2] - min_loc_sorted[0, 2]

    # unique_frames = np.unique(min_loc[:, 2])
    # filtered_points = np.zeros((unique_frames.shape[0], 3))
    # for i, idx in zip(unique_frames, range(unique_frames.shape[0])):
    #     points = min_loc[min_loc[:, 2] == i, :]
    #     centroid = np.round(
    #         np.array([sum(points[:, 0]) / len(points), sum(points[:, 1]) / len(points)])
    #     )
    #     filtered_points[idx, :2] = centroid
    #     filtered_points[idx, 2] = i

    # Find a first frame that is inside the soma (at the minimum)
    filtered_points = min_loc_sorted
    filtered_points = filtered_points[filtered_points[:, 2] >= 0, :]
    points_frame_zero = filtered_points[filtered_points[:, 2] == 0, :]
    points_frame_zero = points_frame_zero[np.argsort(points_frame_zero[:, 3]), :]

    # filtered_points[filtered_points[:, 2] == 0, :] = points_frame_zero

    # Find a last frame that is outside the soma
    _, soma_positions = ident_soma(sta_re, position, radius=soma_radius)
    axon_points = filtered_points[
        ~np.logical_and(
            np.in1d(filtered_points[:, 0], soma_positions[:, 0]),
            np.in1d(filtered_points[:, 1], soma_positions[:, 1]),
        ),
        :,
    ]

    # Build "new" neuron containing one soma position and many axon positions

    points_to_trace = np.zeros((np.shape(axon_points)[0] + 1, 4), dtype=int)
    points_to_trace[0, :] = points_frame_zero[0, :]
    points_to_trace[1:, :] = axon_points

    possible_endpoints = np.where(
        points_to_trace[:, 2] == np.max(points_to_trace[:, 2])
    )[0]

    n = np.max(possible_endpoints) + 1
    pos = {i: (points_to_trace[i, 0], points_to_trace[i, 1]) for i in range(n)}

    length_results = np.empty(len(possible_endpoints), dtype=object)
    shortest_results = np.empty(len(possible_endpoints), dtype=object)
    circles = np.zeros(len(possible_endpoints))

    for endpoint, test in zip(possible_endpoints, range(len(possible_endpoints))):

        circle = 1
        while True:

            try:
                G = nx.random_geometric_graph(n, circle, pos=pos)
                # print(G.nodes)
                shortest_results[test] = nx.shortest_path(G, 0, endpoint)
                length_results[test] = nx.shortest_path_length(G, 0, endpoint)
                circles[test] = circle
                break
            except nx.exception.NetworkXNoPath:
                circle = circle + 1

    # Check if the shortest route is meaningfull

    decision = np.argmin(circles)

    direct_distance = np.sqrt(
        (points_to_trace[0, 0] - points_to_trace[-1, 0]) ** 2
    ) + np.sqrt((points_to_trace[0, 1] - points_to_trace[-1, 1]) ** 2)

    shortest_route = points_to_trace[shortest_results[decision], :]
    length = length_results[decision]
    return filtered_points, shortest_route, length, np.min(circles), direct_distance


def calculate_distance_speed(
    route, sampling_freq=17852.76785, mm_per_channel=2.67 / 64
):

    distances = np.zeros(len(route) - 1)
    for points in range(len(route) - 1):
        distances[points] = np.sqrt(
            (route[points, 0] - route[points + 1, 0]) ** 2
        ) + np.sqrt((route[points, 1] - route[points + 1, 1]) ** 2)
    distance = np.sum(distances) * mm_per_channel
    speed = (distance / (route[-1, 2] * (1 / sampling_freq))) / 1000

    control_distance = np.sqrt((route[0, 0] - route[-1, 0]) ** 2) + np.sqrt(
        (route[0, 1] - route[-1, 1]) ** 2
    )
    if control_distance < 5 * (2.67 / 64):
        speed = 0
    # if len(route) < 3:
    #     speed = 0
    # if route[-1, 2] < 5:
    #     speed = 0

    return distance, speed


def ident_soma(sta_re, positions, radius=4):
    tree = KDTree(positions[:, :2])
    channel = sta_re.shape[0]
    soma = tree.query_ball_point(
        np.unravel_index(np.argmin(np.min(sta_re, axis=2)), (channel, channel)), radius
    )
    soma_size = len(soma)
    return soma_size, positions[soma, :]


class Neuron:
    def __init__(self, sta):
        self.data = {}
        self.data["sta"] = sta.copy()
        self.data["sta_min"] = np.min(sta, axis=1)
        self.tracing = {}

    def imshow(self, dataset="sta_min"):
        fig = px.imshow(self._re(dataset), color_continuous_scale="gray")
        return fig

    def filter_neuron(
        self,
        filter=False,
        filter_array=None,
        threshold=0,
        nr_filter_steps=2,
        neighbours=2,
        filter_threshold=0.2,
        max_channel_nr=400,
    ):
        self.data["positions"], self.data["filtered"] = filter_neuron(
            self.data["sta"],
            filter=filter,
            filter_array=filter_array,
            threshold=threshold,
            nr_filter_steps=nr_filter_steps,
            neighbours=neighbours,
            filter_threshold=filter_threshold,
            max_channel_nr=max_channel_nr,
        )

    def _re(self, dataset="sta"):
        channel = np.sqrt(self.data[dataset].shape[0]).astype(int)

        if self.data[dataset].ndim == 2:
            return np.reshape(
                self.data[dataset], (channel, channel, self.data[dataset].shape[1])
            )
        elif self.data[dataset].ndim == 1:
            return np.reshape(self.data[dataset], (channel, channel))

    def plot_facet(self, dataname="sta"):
        fig = px.imshow(
            self._re(dataname), facet_col=2, facet_col_wrap=5, facet_col_spacing=0.1
        )
        return fig

    def video(self, dataset="sta"):
        fig, ax = plt.subplots(figsize=(2, 2), dpi=250)
        ims = []
        data_re = self._re(dataset=dataset)
        for i in range(self.data[dataset].shape[1]):
            im = ax.imshow(
                data_re[:, :, i],
                animated=True,
                cmap="magma",
                vmin=np.min(self.data[dataset]),
                vmax=np.max(self.data[dataset]),
            )
            # if i == 0:
            #     ax.imshow(
            #         data_re[:, :, 0],
            #         cmap="magma",
            #         vmin=np.min(self.data[dataset]),
            #         vmax=np.max(self.data[dataset]),
            #     )  # show an initial one first
            ims.append([im])
        self.ani = animation.ArtistAnimation(
            fig, ims, interval=50, blit=True, repeat_delay=1000
        )
        return HTML(self.ani.to_jshtml())

    def positions_plot(self):
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=self.data["positions"][:, 1],
                y=self.data["positions"][:, 0],
                mode="markers",
                marker=dict(symbol="square", color="black", colorscale="gray"),
            )
        )
        fig.update_layout(height=500, width=500)
        fig.update_yaxes(range=(64, 0))
        fig.update_xaxes(range=(0, 64))
        return fig

    def interactive_plot(self, window, compare_figure=None):
        interactive_plot = Electric_neuron(
            self._re(), self.positions_plot(), window, compare_figure
        )
        return interactive_plot.return_figure()

    def compare_results_plot(self, window, compare_data="sta_min"):
        min_figure = go.FigureWidget(self.imshow(compare_data))
        scatter_fig = self.interactive_plot(window, compare_figure=min_figure)
        return HBox([scatter_fig, min_figure])

    def trace_neuron(
        self, point_radius=2, bins="auto", use_filter=True, soma_radius=10
    ):
        sta_re = self._re()
        if use_filter:
            sta_new = np.zeros_like(sta_re)
            sta_new[
                self.data["positions"][:, 0], self.data["positions"][:, 1], :
            ] = sta_re[self.data["positions"][:, 0], self.data["positions"][:, 1]]
        else:
            sta_new = sta_re

        self.tracing["filtered_points"], self.tracing["shortest_route"], self.tracing[
            "length"
        ], self.tracing["circle"], self.tracing["direct_distance"] = trace_neuron(
            sta_new,
            self.data["positions"],
            point_radius=point_radius,
            bins=bins,
            soma_radius=soma_radius,
        )

    def filter_and_trace(
        self,
        nr_filter_steps=2,
        neighbours=2,
        filter_threshold=0.2,
        max_channel_nr=400,
        point_radius=2,
    ):
        self.filter_neuron(
            self,
            nr_filter_steps=nr_filter_steps,
            neighbours=neighbours,
            filter_threshold=filter_threshold,
            max_channel_nr=max_channel_nr,
        )
        self.trace_neuron(point_radius=point_radius)
        self.calculate_speed(self)

    def calculate_speed(self):
        self.tracing["distance"], self.tracing["speed"] = calculate_distance_speed(
            self.tracing["shortest_route"]
        )

    def plot_tracing(self, background="filtered"):

        sta_re = self._re()

        sta_new = np.zeros_like(sta_re)
        sta_new[self.data["positions"][:, 0], self.data["positions"][:, 1], :] = sta_re[
            self.data["positions"][:, 0], self.data["positions"][:, 1]
        ]

        fig = self.imshow(background)
        # fig.update_coloraxes(
        #     patch=dict(
        #         colorscale=[
        #             [0, "rgb(250, 250, 250)"],
        #             [0.01, "rgb(200, 200, 200)"],
        #             [0.1, "rgb(150, 150, 150)"],
        #             [0.5, "rgb(100, 100, 100)"],
        #             [0.8, "rgb(50, 50, 50)"],
        #             [1.0, "rgb(0, 0, 0)"],
        #         ]
        #     )
        # )

        for i in range(self.tracing["filtered_points"].shape[0]):
            fig.add_trace(
                go.Scatter(
                    x=[self.tracing["filtered_points"][i, 1]],
                    y=[self.tracing["filtered_points"][i, 0]],
                    mode="markers+text",
                    marker=dict(color="red", opacity=0.2),
                    text=str(self.tracing["filtered_points"][i, 2]),
                    textfont=dict(color="white"),
                )
            )
        fig.add_trace(
            go.Scatter(
                x=self.tracing["shortest_route"][:, 1],
                y=self.tracing["shortest_route"][:, 0],
                mode="lines",
                line=dict(color="red"),
            )
        )

        fig.update_layout(height=1000)
        return fig

    def ident_soma(self, radius=4):
        tree = KDTree(self.data["positions"][:, :2])
        channel = np.sqrt(self.data["sta"].shape[0]).astype(int)
        self.tracing["soma"] = tree.query_ball_point(
            np.unravel_index(np.argmin(np.min(self._re(), axis=2)), (channel, channel)),
            radius,
        )

    def plot_soma(self):
        fig = self.positions_plot()
        fig.add_scatter(
            x=self.data["positions"][self.tracing["soma"], 1],
            y=self.data["positions"][self.tracing["soma"], 0],
            mode="markers",
            marker=dict(color="red"),
        )
        return fig


class Electric_neuron:
    def __init__(self, sta_re, plot, window, compare_figure=None):
        self.reset_figure = True
        self.sta_re = sta_re
        self.window = window
        self.plot = plot
        self.locations_plt = go.FigureWidget(plot)
        self.scatter = self.locations_plt.data[0]
        self.compare_figure = compare_figure

        self.scatter.on_click(self.get_point_channel)

    def get_point_channel(self, trace, points, selector):

        if not selector.ctrl:

            self.channel_indices = []
        # Transforming the string containing the cell index into integer
        channel_ids = np.array([points.xs[0], points.ys[0]])
        # append to the list of cell indices
        self.channel_indices.append(channel_ids)
        color = lambda: random.randint(0, 255)
        hexcolor = "#%02X%02X%02X" % (color(), color(), color())
        self.locations_plt.add_scatter(
            x=points.xs,
            y=points.ys,
            mode="markers",
            showlegend=False,
            hoverinfo="skip",
            marker=dict(color=hexcolor, size=5),
        )

        if type(self.compare_figure) == type(go.FigureWidget()):
            print("yes")
            self.compare_figure.add_scatter(
                x=points.xs,
                y=points.ys,
                mode="markers",
                showlegend=False,
                hoverinfo="skip",
                marker=dict(color=hexcolor, size=5),
            )
        # create scatterplot to highlight selected points

        # clears previous output
        self.window.clear_output()
        # Call function that provides waveforms
        self.plot_channels(hexcolor)

    def return_figure(self):
        return self.locations_plt

    def plot_channels(self, hexcolor):
        if self.reset_figure:
            self.fig = go.Figure()

            self.reset_figure = False
        for channel in range(len(self.channel_indices)):
            name = (
                str(self.channel_indices[channel][0])
                + " "
                + str(self.channel_indices[channel][1])
            )

            self.fig.add_trace(
                go.Scatter(
                    y=self.sta_re[
                        self.channel_indices[channel][1],
                        self.channel_indices[channel][0],
                        :,
                    ],
                    name=name,
                    mode="lines",
                    line=dict(color=hexcolor),
                )
            )
        with self.window:
            display(self.fig)
