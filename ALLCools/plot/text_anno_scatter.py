import pandas as pd
import matplotlib
import numpy as np

def _calculate_luminance(color):
	"""
	Calculate the relative luminance of a color according to W3C standards

	Parameters
	----------
	color : matplotlib color or sequence of matplotlib colors
		Hex code, rgb-tuple, or html color name.
	Returns
	-------
	luminance : float(s) between 0 and 1

	"""
	rgb = matplotlib.colors.colorConverter.to_rgba_array(color)[:, :3]
	rgb = np.where(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
	lum = rgb.dot([0.2126, 0.7152, 0.0722])
	try:
		return lum.item()
	except ValueError:
		return lum


def _text_anno_scatter(
    data: pd.DataFrame,
    ax,
    x: str,
    y: str,
    palette: dict = None,
    dodge_text=False,
    anno_col="text_anno",
    text_anno_kws=None,
    text_transform=None,
    dodge_kws=None,
    linewidth=0.5,
    labelsize=5,
    bbox_alpha=1,
    luminance=0.48
):
    """Add text annotation to a scatter plot."""
    # prepare kws
    _text_anno_kws = {
        "fontsize": labelsize,
        "fontweight": "black",
        "horizontalalignment": "center",
        "verticalalignment": "center",
        'edge_color':(0.5, 0.5, 0.5, 0.2),
        'face_color':(0.8, 0.8, 0.8, 0.2),
    }
    if text_anno_kws is None:
        text_anno_kws={}
    _text_anno_kws.update(text_anno_kws)
    face_color=_text_anno_kws.pop('face_color')
    edge_color=_text_anno_kws.pop('edge_color')
    
    colors=text_anno_kws.pop('color','black')
    # plot each text
    text_list = []
    for text, sub_df in data.groupby(anno_col):
        if text_transform is None:
            text = str(text)
        else:
            text = text_transform(text)
        if text.lower() in ["", "nan"]:
            continue
        _x, _y = sub_df[[x, y]].median()

        if palette is not None:
            _fc = palette[text]
        else:
            _fc = face_color
        if type(colors)==dict:
            color=colors[text]
        else:
            color=colors
        lum = _calculate_luminance(color)
        if lum > luminance:
            color='black'
        if _fc is None:
            bbox=None
        else:
            bbox={"boxstyle": "round", "ec": edge_color, "fc": _fc, "linewidth": linewidth,'alpha':bbox_alpha}
        text = ax.text(
            _x,
            _y,
            text,
            color=color,
            fontdict=_text_anno_kws,
            bbox=bbox,
        )
        text_list.append(text)

    if dodge_text:
        try:
            from adjustText import adjust_text

            _dodge_parms = {
                "force_points": (0.02, 0.05),
                "arrowprops": {
                    "arrowstyle": "->",
                    "fc": edge_color,
                    "ec": "none",
                    "connectionstyle": "angle,angleA=-90,angleB=180,rad=5",
                },
                "autoalign": "xy",
            }
            if dodge_kws is not None:
                _dodge_parms.update(dodge_kws)
            adjust_text(text_list, x=data["x"], y=data["y"], **_dodge_parms)
        except ModuleNotFoundError:
            print("Install adjustText package to dodge text, see its github page for help")

    return
