function [c1,c2,pc_DMN,pc_TPN,DMN_f,yp_DMN,TPN_f,yp_TPN,yc_DMN,yc_TPN] ...
    = iter_optimal_PSTP(ITP1,DMN_f_all,...
          yp_DMN_all,pc_DMN_all,pc_TPN_all,TPN_f_all,yp_TPN_all,...
          yc_DMN_all,yc_TPN_all,c1_all,c2_all)


DMN_f = DMN_f_all{ITP1};
yp_DMN = yp_DMN_all(ITP1,:);

pc_DMN = pc_DMN_all(ITP1);
pc_TPN = pc_TPN_all(ITP1);

TPN_f = TPN_f_all{ITP1};
yp_TPN = yp_TPN_all(ITP1,:);

yc_DMN = yc_DMN_all(ITP1,:);
yc_TPN = yc_TPN_all(ITP1,:);

c1 = c1_all(ITP1);
c2 = c2_all(ITP1);


end

