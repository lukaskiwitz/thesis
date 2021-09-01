import os
import sys
from functools import reduce

import numpy as np
from numpy import sin, cos
from paraview.simple import *


def load_cytokine_xdmf(file_path):
    # create a new 'Xdmf3ReaderS'
    field_xdmf = Xdmf3ReaderS(FileName=[file_path])

    return field_xdmf


def load_markers_xdmf(file_path):
    markers_0xdmf = Xdmf3ReaderS(FileName=[file_path])

    return markers_0xdmf


def marker_threshold(source, view, range=(1, 1)):
    threshold = Threshold(Input=source)
    threshold.Scalars = ['CELLS', source.CellData.GetArray(0).Name]
    threshold.ThresholdRange = range

    return threshold


def create_slice(source, origin=(0, 0, 0), normal=(0, 0, 1)):
    slice = Slice(Input=source)
    slice.SliceType = 'Plane'
    slice.HyperTreeGridSlicer = 'Plane'
    slice.SliceOffsetValues = [0.0]
    slice.SliceType.Origin = origin
    slice.HyperTreeGridSlicer.Origin = origin
    slice.SliceType.Normal = normal

    return slice


def create_calc(source, factor="1e18", result_name="Cytokine"):
    calc = Calculator(Input=source)

    calc.Function = '{name}*{f}'.format(f=factor, name=calc.PointData.GetArray(0).Name)
    calc.ResultArrayName = result_name

    Hide(calc)
    return calc


def save_render_view(render_view, file_path, settings):
    render_view.OrientationAxesVisibility = settings["orientation_axes"]
    LoadPalette(paletteName='WhiteBackground')
    render_view.Update()

    SaveScreenshot(file_path, render_view)


def export_render_view(render_view, file_path, settings):
    render_view.OrientationAxesVisibility = settings["orientation_axes"]
    LoadPalette(paletteName='WhiteBackground')
    render_view.Update()

    ExportView(file_path, view=render_view)


def format_cytokine_bar(arrayname, view, settings):
    LUT = GetColorTransferFunction(arrayname)
    PWF = GetOpacityTransferFunction(arrayname)

    LUT.ApplyPreset(settings["field_color_preset"], True)

    color_range = settings["color_bar_range"]

    LUT.RescaleTransferFunction(color_range[0], color_range[1])
    if "opacity_range" in settings:
        opacity_range = settings["opacity_range"]
        PWF.RescaleTransferFunction(opacity_range[0], opacity_range[1])

    LUTColorBar = GetScalarBar(LUT, view)
    LUTColorBar.AutoOrient = 0

    LUTColorBar.WindowLocation = 'AnyLocation'
    LUTColorBar.Position = [0.5, 0.1]
    LUTColorBar.DrawTickMarks = 1
    LUTColorBar.DrawTickLabels = 1
    LUTColorBar.AutomaticLabelFormat = 0

    LUTColorBar.UseCustomLabels = 1
    LUTColorBar.CustomLabels = np.linspace(color_range[0], color_range[1], 5)

    LUTColorBar.LabelFormat = settings["number_format"]
    LUTColorBar.RangeLabelFormat = settings["number_format"]


def format_marker_bar(arrayname, view, cell_type_lookup):
    annotations = [[str(k), str(v[0])] for k, v in cell_type_lookup.items()]
    indexed_colors = [[v[1][0], v[1][1], v[1][2]] for k, v in cell_type_lookup.items()]
    indexed_colors = reduce(lambda x, y: x + y, indexed_colors)
    annotations = reduce(lambda x, y: x + y, annotations)

    fLUT = GetColorTransferFunction(arrayname)

    fLUT.InterpretValuesAsCategories = 1
    fLUT.AnnotationsInitialized = 0
    fLUT.Annotations = annotations
    fLUT.IndexedColors = indexed_colors
    # fLUT.IndexedOpacities = [1.0] * len(cell_type_lookup)

    fLUTColorBar = GetScalarBar(fLUT, view)
    fLUTColorBar.AutoOrient = 0
    fLUTColorBar.WindowLocation = 'AnyLocation'
    fLUTColorBar.Position = [0.5, 0.1]

    fLUTColorBar.DrawAnnotations = 1

    fLUTColorBar.Title = "Cell Type"
    fLUTColorBar.TitleFontSize = settings["cell_type_title_font_size"]
    fLUTColorBar.LabelFontSize = settings["cell_type_label_font_size"]


def render_slice(slice_view, slice_display, origin, img_path, settings):
    slice_display.DataAxesGrid.GridAxesVisibility = 1
    slice_display.DataAxesGrid.XTitle = 'X $(\mu m)$'
    slice_display.DataAxesGrid.YTitle = ''
    slice_display.DataAxesGrid.XTitleFontSize = settings["axis_title_font_size"]
    slice_display.DataAxesGrid.YTitleFontSize = settings["axis_title_font_size"]
    slice_display.DataAxesGrid.XLabelFontSize = settings["axis_label_font_size"]
    slice_display.DataAxesGrid.YLabelFontSize = settings["axis_label_font_size"]
    slice_display.DataAxesGrid.XAxisNotation = 'Fixed'
    slice_display.DataAxesGrid.XAxisPrecision = 0
    slice_display.DataAxesGrid.YAxisNotation = 'Fixed'
    slice_display.DataAxesGrid.YAxisPrecision = 0
    slice_display.DataAxesGrid.AxesToLabel = 3
    slice_display.DataAxesGrid.ShowTicks = settings["axis_ticks"]
    slice_display.DataAxesGrid.ShowEdges = settings["axis_edges"]

    slice_view.EnableRayTracing = 0
    slice_view.InteractionMode = '2D'
    slice_view.CameraPosition = origin + np.array(settings["slice_normal"])
    slice_view.CameraFocalPoint = origin
    slice_view.CameraViewUp = [0, 1, 0]
    slice_view.CameraParallelScale = 100 * settings["slice_zoom"]
    slice_view.ViewSize = settings["render_view_size"]

    export_render_view(slice_view, os.path.join(img_path, "slice.pdf"), settings)
    save_render_view(slice_view, os.path.join(img_path, "slice.png"), settings)


def render_volume(volume_view, origin, img_path, settings):
    volume_view.EnableRayTracing = settings["volume_raytracing"]
    volume_view.ProgressivePasses = settings["volume_raytracing_progressive_passes"]
    volume_view.SamplesPerPixel = settings["volume_raytracing_samples"]
    volume_view.AmbientSamples = settings["volume_raytracing_ambient_samples"]
    volume_view.LightScale = settings["volume_raytracing_light_scale"]
    volume_view.InteractionMode = '3D'
    volume_view.CameraPosition = settings["volume_camera_pos"]
    volume_view.CameraFocalPoint = origin
    volume_view.CameraViewUp = [0, 0, 1]

    volume_view.ViewSize = settings["render_view_size"]
    save_render_view(volume_view, os.path.join(img_path, "volume.png"), settings)


def render_markers(marker_view, origin, img_path, settings):
    marker_view.EnableRayTracing = settings["volume_raytracing"]
    marker_view.ProgressivePasses = settings["volume_raytracing_progressive_passes"]
    marker_view.SamplesPerPixel = settings["volume_raytracing_samples"]
    marker_view.InteractionMode = '3D'
    marker_view.CameraPosition = settings["volume_camera_pos"]
    marker_view.CameraFocalPoint = origin
    marker_view.CameraViewUp = [0, 0, 1]

    marker_view.ViewSize = settings["render_view_size"]
    save_render_view(marker_view, os.path.join(img_path, "marker.png"), settings)


def render_legend(legend_view, img_path, settings, suffix):
    legend_view.ViewSize = settings["render_view_size"]
    export_render_view(legend_view, os.path.join(img_path, "legend_{s}.pdf".format(s=suffix)), settings)


def make_images(cytokine_path, marker_path, img_path, settings):
    volume_view = GetActiveViewOrCreate('RenderView')  # render volume_view
    slice_view = CreateView('RenderView')  # slice volume_view
    legend_view_1 = CreateView('RenderView')  # legend volume_view
    legend_view_2 = CreateView('RenderView')  # legend volume_view

    marker_view = CreateView('RenderView')  # only cell markers

    source = load_cytokine_xdmf(cytokine_path)  # loads file
    markers = load_markers_xdmf(marker_path)  # loads boundary markers

    original_field_name = source.PointData.GetArray(0).Name

    source = create_calc(source,
                         factor=settings["conversion_factor"],
                         result_name=settings["field_name"]
                         )  # unit conversion

    bb = source.GetDataInformation().GetBounds()
    origin = np.array(
        [
            np.mean([bb[0], bb[1]]),
            np.mean([bb[2], bb[3]]),
            np.mean([bb[4], bb[5]])
        ])

    r, t, p = settings["volume_camera_pos"]

    x = r * sin(np.deg2rad(t)) * cos(np.deg2rad(p))
    y = r * sin(np.deg2rad(t)) * sin(np.deg2rad(p))
    z = r * cos(np.deg2rad(t))

    settings["volume_camera_pos"] = np.array([x, y, z]) + origin

    #######
    volume_display = GetDisplayProperties(source, view=volume_view)
    volume_display.SetRepresentationType('Volume')
    volume_display.SelectMapper = 'Resample To Image'
    ####################

    cell_type_thresholds = []
    marker_display = None

    for k, v in settings["lookup"].items():
        threshold = marker_threshold(markers, volume_view, range=[int(k), int(k)])
        cell_type_thresholds.append(threshold)
        marker_display = Show(threshold, volume_view, 'UnstructuredGridRepresentation')
        marker_display2 = Show(threshold, marker_view, 'UnstructuredGridRepresentation')

        marker_display2.Opacity = v[2] if settings["volume_raytracing"] == 1 and settings[
            "marker_view_uniform_opacity"] == 0 else 1
        marker_display.Opacity = v[2] if settings["volume_raytracing"] == 1 else 1

    slice = create_slice(source, origin=origin + settings["slice_origin"],
                         normal=settings["slice_normal"])  # creates slice
    slice_display = Show(slice, slice_view, 'GeometryRepresentation')  # slice display

    for threshold in cell_type_thresholds:
        if threshold.CellData.GetArray(0) is None:
            continue  # no cells present

        clip = Clip(Input=threshold)
        clip.ClipType = 'Box'
        clip.HyperTreeGridClipper = 'Plane'
        clip.Scalars = ['CELLS', threshold.CellData.GetArray(0).Name]

        l = np.array([abs(bb[0] - bb[1]), abs(bb[2] - bb[3]), settings["layer_distance"]])

        clip.ClipType.Position = np.array(origin + settings["slice_origin"]) - 0.5 * l - np.array(
            settings["slice_normal"]) * settings["layer_distance"] / 2
        clip.ClipType.Length = l

        display = Show(clip, slice_view, 'UnstructuredGridRepresentation')
        display.Opacity = 1

    volume_display.SetScalarBarVisibility(legend_view_1, True)

    if marker_display is not None:
        marker_display.SetScalarBarVisibility(legend_view_2, True)

    output_array = source.PointData.GetArray(0)

    for i in range(source.PointData.GetNumberOfArrays()):
        if source.PointData.GetArray(i).Name == settings["field_name"]:
            output_array = source.PointData.GetArray(i)

    if "color_bar_range" not in settings:
        settings["color_bar_range"] = output_array.GetRange()

    elif len(settings["color_bar_range"]) == 1:
        settings["color_bar_range"] = [settings["color_bar_range"][0], output_array.GetRange()[1]]

    format_cytokine_bar(output_array.Name, legend_view_1, settings)

    format_marker_bar(markers.CellData.GetArray(0).Name, legend_view_2, settings["lookup"])

    render_slice(slice_view, slice_display, origin, img_path, settings)
    render_volume(volume_view, origin, img_path, settings)
    render_markers(marker_view, origin, img_path, settings)
    render_legend(legend_view_1, img_path, settings, "cytokine")
    render_legend(legend_view_2, img_path, settings, "marker")


settings = {
    "orientation_axes": False,
    "render_view_size": [400, 400],
    "conversion_factor": 1e15,
    "field_name": "Cytokine(nM)",
    "cell_type_title_font_size": 20,
    "cell_type_label_font_size": 20,
    "slice_zoom": 1.2,
    "slice_origin": [0, 0, 100],
    "slice_normal": [0, 0, 1],
    "layer_distance": 20,
    "axis_title_font_size": 15,
    "axis_label_font_size": 15,
    "number_format": '%2.1g',
    "axis_ticks": 1,
    "axis_edges": 0,
    "field_color_preset": "Cool to Warm",
    "volume_camera_pos": [850, 70, -110],
    "volume_raytracing": 0,
    "volume_raytracing_progressive_passes": 0,
    "volume_raytracing_samples": 2,
    "volume_raytracing_ambient_samples": 2,
    "volume_raytracing_light_scale": 1,
    "marker_view_uniform_opacity": 1,
    "lookup": {
        "1": ["default", (0.5, 0.5, 0.5), 1],
    }
}

if __name__ == "__main__":
    import json

    xmdf = sys.argv[1]
    marker_path = sys.argv[2]
    img_path = sys.argv[3]

    if len(sys.argv) > 4:
        settings_path = sys.argv[4]
        with open(settings_path, "r") as f:
            settings.update(json.load(f))

    make_images(xmdf, marker_path, img_path, settings)
