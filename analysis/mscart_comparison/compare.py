"""Plot MUSS-MSCART comparisons from the saved ground and satellite data."""
from pathlib import Path
import argparse
import csv
import json
import os
import numpy as np
from scipy.io import netcdf_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

ROOT = Path(__file__).resolve().parent


def csv_write(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_model(path, level, fov, depth):
    with path.open(encoding='utf-8-sig', newline='') as f:
        rows = [r for r in csv.DictReader(f)
                if r['extinction'] == level and float(r['fov_full_mrad']) == fov]
    channels, orders = [], []
    for repeat in sorted({r['repeat'] for r in rows}):
        rr = sorted([r for r in rows if r['repeat'] == repeat], key=lambda r: float(r['depth_m']))
        np.testing.assert_allclose([float(r['depth_m']) for r in rr], depth, atol=1e-10, rtol=0)
        channels.append([[float(r['parallel']), float(r['perpendicular'])] for r in rr])
        orders.append([[[float(r[f'order{k}_{c}']) for c in ['parallel', 'perpendicular']]
                         for k in [1, 2, 3, 4]] for r in rr])
    assert channels, f'Missing model case: {level}, {fov}'
    channels, orders = np.asarray(channels), np.asarray(orders)
    assert np.isfinite(channels).all() and np.isfinite(orders).all()
    np.testing.assert_allclose(orders.sum(2), channels, rtol=1e-10, atol=1e-21)
    return channels, orders


def read_mscart(path, fov, ground):
    with netcdf_file(path, 'r', mmap=False) as ds:
        fovs = ds.variables['fov'][:].astype(float)
        index = int(np.argmin(abs(fovs*1000-fov)))
        np.testing.assert_allclose(fovs[index]*1000, fov, rtol=1e-6)
        total = ds.variables['RunMeanRCRad'][:].astype(float)[:, 0, index]
        orders = ds.variables['RunMeanNthRCRad'][:].astype(float)[:, 0, index]
        np.testing.assert_allclose(total.mean(0), ds.variables['MeanRCRad'][:][0, index], rtol=1e-14, atol=0)
    if ground:
        total /= 1e6
        orders /= 1e6
    channels = np.stack([(total[..., 0]+total[..., 1])/2, (total[..., 0]-total[..., 1])/2], axis=-1)
    nth = np.stack([(orders[..., 0]+orders[..., 1])/2, (orders[..., 0]-orders[..., 1])/2], axis=-1)
    nth = np.concatenate([nth, (channels-nth.sum(2))[:, :, None]], axis=2)
    np.testing.assert_allclose(nth.sum(2), channels, rtol=1e-10, atol=1e-21)
    return channels, nth


class SingleOriginZero(ScalarFormatter):
    def __call__(self, value, pos=None):
        return '' if value == 0 else super().__call__(value, pos)


def style():
    font = Path('C:/Windows/Fonts/msyh.ttc')
    if font.exists(): fontManager.addfont(str(font))
    family = FontProperties(fname=str(font)).get_name() if font.exists() else ['Noto Sans CJK SC','WenQuanYi Zen Hei','DejaVu Sans']
    plt.rcParams.update({'font.family': family, 'font.size': 9, 'axes.labelsize': 9,
        'axes.titlesize': 9, 'legend.fontsize': 9, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
        'axes.linewidth': .65, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'mathtext.fontset': 'dejavusans', 'pdf.fonttype': 42, 'axes.unicode_minus': False})


def divide(n, d):
    return np.divide(n, d, out=np.full_like(n, np.nan), where=d > 0)


def compare_case(case, model_file=None):
    ground = case == 'ground'
    folder = ROOT / case
    if model_file is None:
        model_file = folder/'model_results.csv'
    model_file = model_file.resolve()
    print(f'[{case}] Reading model CSV: {model_file}', flush=True)
    depth = (np.arange(60 if ground else 100)+.5)*(5 if ground else 3)
    edge = depth+(2.5 if ground else 1.5)
    configs = [(level, fov) for level in ['low', 'high'] for fov in [.5, 2.]] if ground else [('hu2001', fov) for fov in [.001, .13, 13.]]
    style()
    fig, axes = plt.subplots(2 if ground else 1, 2 if ground else 3,
                             figsize=(14/2.54, (12.4 if ground else 7.5)/2.54), squeeze=False)
    fig.subplots_adjust(left=.16 if ground else .13, right=.965, bottom=.165 if ground else .26,
                        top=.79 if ground else .70, wspace=.38, hspace=.80)
    fig.suptitle('本模型与 MSCART：'+('地基案例' if ground else '星载案例'), fontsize=10, y=.98)
    handles = [Line2D([], [], color='#0072B2', lw=1.25),
               Line2D([], [], color='#222222', lw=.95, ls='--', marker='o', ms=3, mfc='white')]
    fig.legend(handles, ['本模型', 'MSCART'], loc='upper center', bbox_to_anchor=(.54,.925),
               ncol=2, frameon=False, handlelength=2.1, columnspacing=2.4)
    rows, summaries = [], []
    for panel, (level, fov) in enumerate(configs):
        model, _ = read_model(model_file, level, fov, depth)
        mscart, _ = read_mscart(folder/(f'{level}_result.nc' if ground else 'result.nc'), fov, ground)
        ax = axes.flat[panel]
        display = edge <= (150 if ground else 300)
        for label, channels, color, marker in [('model',model,'#0072B2',False),('MSCART',mscart,'#222222',True)]:
            mean = channels.mean(0)
            integrated = divide(np.cumsum(mean[:,1]), np.cumsum(mean[:,0]))
            repeats = divide(np.cumsum(channels[:,:,1],axis=1),np.cumsum(channels[:,:,0],axis=1))
            ax.fill_between(edge[display], repeats.min(0)[display], repeats.max(0)[display],color=color,alpha=.11,lw=0)
            ax.plot(edge[display], integrated[display],color=color,lw=.95 if marker else 1.25,
                    ls='--' if marker else '-',marker='o' if marker else None,ms=2.8,
                    mfc='white',mew=.7,markevery=6 if ground else 10)
            for j in range(len(depth)):
                rows.append(dict(extinction=level,fov_full_mrad=fov,solver=label,depth_m=float(depth[j]),
                    integration_end_m=float(edge[j]),parallel=float(mean[j,0]),perpendicular=float(mean[j,1]),
                    delta_pointwise=float(divide(mean[:,1],mean[:,0])[j]),delta_cumulative=float(integrated[j])))
        windows = [(20,120)] if ground else [(0,150),(0,300)]
        for start,end in windows:
            mask=(depth>=start)&(depth<end)
            m=model.mean(0)[mask].sum(0);r=mscart.mean(0)[mask].sum(0)
            dm=m[1]/m[0];dr=r[1]/r[0]
            summaries.append(dict(extinction=level,fov_full_mrad=fov,start_m=start,end_m=end,
                model_delta=float(dm),mscart_delta=float(dr),relative_difference_percent=float(100*(dm/dr-1)),
                parallel_difference_percent=float(100*(m[0]/r[0]-1)),perpendicular_difference_percent=float(100*(m[1]/r[1]-1))))
        ax.set_title(('低消光' if level=='low' else '高消光' if level=='high' else '全视场')+f' / {fov:g} mrad',pad=6)
        ax.set_xlim(0,150 if ground else 300)
        ax.set_xticks([0,50,100,150] if ground else [0,100,200,300])
        ax.set_ylim(bottom=0)
        ax.set_xlabel('云内积分距离 (m)')
        if ground or panel==0: ax.set_ylabel(r'累计线退偏比 $\delta(D)$')
        fmt=SingleOriginZero(useMathText=True);fmt.set_scientific(False)
        ax.yaxis.set_major_formatter(fmt)
        ax.text(.04,.95,f'({chr(97+panel)})',transform=ax.transAxes,fontweight='bold',va='top',zorder=10,
                bbox=dict(facecolor='white',edgecolor='none',pad=.2))
        ax.grid(color='#e7e7e7',lw=.4);ax.set_axisbelow(True)
    fig.text(.16 if ground else .13,.035 if ground else .065,'通道先平均，再积分取比；阴影为独立重复范围。',fontsize=8)
    for suffix in ['png','pdf']:
        fig.savefig(folder/f'comparison.{suffix}',dpi=600,facecolor='white')
    plt.close(fig)
    csv_write(folder/'comparison.csv',rows)
    csv_write(folder/'summary.csv',summaries)
    image_file = (folder/'comparison.png').resolve()
    print(f'[{case}] Saved PNG: {image_file}', flush=True)
    print(f'[{case}] Saved PDF: {(folder/"comparison.pdf").resolve()}', flush=True)
    print(f'[{case}] Saved tables: comparison.csv, summary.csv', flush=True)
    print(json.dumps(dict(case=case,configurations=len(configs),source_model=str(model_file),
                         summary_rows=len(summaries)),ensure_ascii=False), flush=True)
    return image_file


def open_generated_image(image_file):
    if os.name != 'nt':
        print(f'Open the saved image to view the plot: {image_file}', flush=True)
        return
    try:
        os.startfile(str(image_file))
        print(f'Opened image viewer: {image_file}', flush=True)
    except OSError as error:
        print(f'Image was saved; could not open the viewer: {error}\n{image_file}', flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=['ground','satellite','all'],default='all')
    p.add_argument('--model-file',type=Path,help='Path to a saved MUSS result CSV; use with a single --case.')
    p.add_argument('--no-show',action='store_true',help='Save output files without opening the image viewer.')
    args=p.parse_args()
    if args.model_file and args.case=='all': p.error('--model-file requires --case ground or satellite')
    if args.model_file:
        candidate = args.model_file.expanduser()
        if not candidate.is_absolute() and not candidate.is_file():
            candidate = ROOT/candidate
        args.model_file = candidate.resolve()
    generated = []
    for case in ['ground','satellite'] if args.case=='all' else [args.case]:
        generated.append(compare_case(case,args.model_file))
    if not args.no_show:
        for image_file in generated:
            open_generated_image(image_file)


if __name__=='__main__': main()
