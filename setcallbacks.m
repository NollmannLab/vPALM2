%% vPALM set callbacks    
function setcallbacks(h) 

set(h.push_load_file,'callback',{@vPALM_load_file,h});
set(h.plot_localizations_3d,'callback',{@vPALM_plot_localizations_3d,h});
set(h.plot_localizations_2d,'callback',{@vPALM_plot_localizations_2d,h});
set(h.plot_reconstruction,'callback',{@vPALM_plot_reconstruction,h});
set(h.plot_density,'callback',{@plot_density,h});
set(h.cal_3d,'callback',{@vPALM_cal_3d,h});
set(h.calc_res,'callback',{@vPALM_calc_res,h});
set(h.Zcolormap,'callback',{@Zcolormap,h});
set(h.zplotoption,'callback',{@zplotoption,h});
set(h.stats,'callback',{@vPALM_stats,h});
set(h.save_data,'callback',{@vPALM_save_data,h});
set(h.apply_3dcal,'callback',{@vPALM_apply_3dcal,h});
set(h.save_PALMvis,'callback',{@vPALM_save_PALMvis,h});
set(h.make_roi,'callback',{@vPALM_make_roi,h});
set(h.reset_roi,'callback',{@vPALM_reset_roi,h});
set(h.save_axes1,'callback',{@vPALM_save_axes1,h});
set(h.offsetx_slider,'callback',{@vPALM_offsetx_slider,h});
set(h.offsety_slider,'callback',{@vPALM_offsety_slider,h});
set(h.selectbeads,'callback',{@vPALM_selectbeads,h});

% set(h.showbeads,'callback',{@showbeads,h});
