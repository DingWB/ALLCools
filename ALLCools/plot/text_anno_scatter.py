import pandas as pd
import matplotlib
import numpy as np
import copy

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
    dodge_text=False,
    anno_col="text_anno",
    text_kws=None,
    text_transform=None,
    dodge_kws=None,
    luminance=0.48
):
    """Add text annotation to a scatter plot."""
    # prepare kws
    text_kws={} if text_kws is None else text_kws
    text_kws.setdefault("fontsize",5)
    text_kws.setdefault("fontweight","black")
    text_kws.setdefault("ha","center") #horizontalalignment
    text_kws.setdefault("va","center") #verticalalignment
    text_kws.setdefault("color","black") #c
    bbox=dict(boxstyle='round',edgecolor=(0.5, 0.5, 0.5, 0.2),fill=False,
                                facecolor=(0.8, 0.8, 0.8, 0.2),alpha=1,linewidth=0.5)
    text_kws.setdefault("bbox",bbox)
    for key in bbox:
        if key not in text_kws['bbox']:
            text_kws['bbox'][key]=bbox[key]
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
        
        use_text_kws=copy.deepcopy(text_kws) #text_kws.copy()
        if isinstance(text_kws['bbox']['facecolor'],dict):
            use_text_kws['bbox']['facecolor']=text_kws['bbox']['facecolor'][text]
        if isinstance(text_kws['color'],dict):
            use_color=text_kws['color'][text]
            use_text_kws['color']=use_color
        if not luminance is None:
            lum = _calculate_luminance(use_text_kws['bbox']['facecolor'])
            if lum > luminance:
                use_text_kws['color']='black'
                use_text_kws['bbox']['edgecolor']='black'
        
        text = ax.text(
            _x,
            _y,
            text,
            **use_text_kws
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
