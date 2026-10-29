function v = fl_project_box_controls(u,~,args)
% Project the floating nodal control values onto the pointwise bounds.
% This projection belongs to the nodal search model.

bounds =  args{14};
v = min(max(bounds),max(min(bounds),u));
end

