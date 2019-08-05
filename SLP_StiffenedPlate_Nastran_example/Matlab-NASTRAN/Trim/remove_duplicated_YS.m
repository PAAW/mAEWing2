
function [YS_updated,ZS_updated,XS_updated,REFC_updated,REFS_updated,RigidLift_updated,ElasticLift_updated,ElasticMoment_updated] = remove_duplicated_YS(YS,ZS,XS,REFC,REFS,RigidLift,ElasticLift,ElasticMoment)


winglet_Y = max(YS);


% 1 - remove Y for winglet
% 2 - find out distance in strip Y in YS_wing less than 0.5
% 3 - combine them for new lift coefficient

winglet_id = find(abs(YS) == winglet_Y);
YS_winglet = YS(winglet_id,:);
ZS_winglet = ZS(winglet_id,:);
XS_winglet = XS(winglet_id,:);
REFC_winglet = REFC(winglet_id,:);
REFS_winglet = REFS(winglet_id,:);
RigidLift_winglet = RigidLift(winglet_id,:);
ElasticLift_winglet = ElasticLift(winglet_id,:);
ElasticMoment_winglet = ElasticMoment(winglet_id,:);

wing_id = setdiff(1:length(YS),winglet_id);

YS_wing = YS(wing_id,:);
ZS_wing = ZS(wing_id,:);
XS_wing = XS(wing_id,:);
REFC_wing = REFC(wing_id,:);
REFS_wing = REFS(wing_id,:);
RigidLift_wing = RigidLift(wing_id,:);
ElasticLift_wing = ElasticLift(wing_id,:);
ElasticMoment_wing = ElasticMoment(wing_id,:);

YS_tol = 0.5;


for ii = 1:length(YS_wing)

    for jj = ii+1:length(YS_wing)-1
         
        diff = YS_wing(jj) - YS_wing(ii);
        ac_diff = ZS_wing(jj) - ZS_wing(ii);
        
        if abs(diff)< YS_tol 
            
            
            area_temp =   (REFS_wing(ii) +REFS_wing(jj) );
            chord_temp =  REFC_wing(ii) + REFC_wing(jj);
            zs_temp = ZS_wing(ii);
            ys_temp = YS_wing(ii);
            xs_temp = XS_wing(ii) - REFC_wing(ii)/4 + chord_temp/4;
            rigid_cn_temp =(RigidLift_wing(ii)*REFS_wing(ii) + RigidLift_wing(jj)*REFS_wing(jj) )/area_temp;
            
            elastic_cn_temp =(ElasticLift_wing(ii)*REFS_wing(ii) + ElasticLift_wing(jj)*REFS_wing(jj) )/area_temp;

            elastic_m_temp =(ElasticMoment_wing(ii)*REFS_wing(ii)* REFC_wing(ii)+ ElasticMoment_wing(jj)*REFS_wing(jj)*REFC_wing(jj))/area_temp/chord_temp;

            % update YS_wing, and others_wing
            
            updated_index = setdiff(1:length(YS_wing),jj);
            
            YS_wing = YS_wing(updated_index,:); YS_wing(ii) = ys_temp;
            ZS_wing = ZS_wing(updated_index,:); ZS_wing(ii) = zs_temp;
            XS_wing = XS_wing(updated_index,:); XS_wing(ii) = xs_temp;
            
            REFC_wing = REFC_wing(updated_index,:); REFC_wing(ii) = chord_temp;
            REFS_wing = REFS_wing(updated_index,:); REFS_wing(ii) = area_temp;
            
            RigidLift_wing = RigidLift_wing(updated_index,:);RigidLift_wing(ii) = rigid_cn_temp;
            ElasticLift_wing = ElasticLift_wing(updated_index,:);ElasticLift_wing(ii) = elastic_cn_temp;
            
            ElasticMoment_wing = ElasticMoment_wing(updated_index,:); ElasticMoment_wing(ii) = elastic_m_temp;
            
        end
    end
end

%%
% id_sequence = [1:length(YS_wing)/2 length(YS_wing)/2+1:length(YS_wing)+length(winglet_Y)/2 ...
%     1:length(YS_wing)/2 length(YS_wing)/2+1:length(YS_wing)+length(winglet_Y)/2]
% 
% 
% updated_strips = length(YS_wing) + length(winglet_Y);

%
YS_updated = [YS_wing(1:length(YS_wing)/2) ; YS_winglet(1:length(winglet_id)/2) ;YS_wing(1+length(YS_wing)/2:length(YS_wing))  ; YS_winglet(1+length(winglet_id)/2:end)];
ZS_updated = [ZS_wing(1:length(ZS_wing)/2) ; ZS_winglet(1:length(winglet_id)/2) ;ZS_wing(1+length(ZS_wing)/2:length(ZS_wing))  ; ZS_winglet(1+length(winglet_id)/2:end)];
XS_updated = [XS_wing(1:length(XS_wing)/2) ; XS_winglet(1:length(winglet_id)/2) ;XS_wing(1+length(XS_wing)/2:length(XS_wing))  ; XS_winglet(1+length(winglet_id)/2:end)];
REFC_updated = [REFC_wing(1:length(REFC_wing)/2) ; REFC_winglet(1:length(winglet_id)/2) ;REFC_wing(1+length(REFC_wing)/2:length(REFC_wing))  ; REFC_winglet(1+length(winglet_id)/2:end)];
REFS_updated = [REFS_wing(1:length(REFS_wing)/2) ; REFS_winglet(1:length(winglet_id)/2) ;REFS_wing(1+length(REFS_wing)/2:length(REFS_wing))  ; REFS_winglet(1+length(winglet_id)/2:end)];
RigidLift_updated = [RigidLift_wing(1:length(RigidLift_wing)/2) ; RigidLift_winglet(1:length(winglet_id)/2) ;RigidLift_wing(1+length(RigidLift_wing)/2:length(RigidLift_wing))  ; RigidLift_winglet(1+length(winglet_id)/2:end)];
ElasticLift_updated = [ElasticLift_wing(1:length(ElasticLift_wing)/2) ; ElasticLift_winglet(1:length(winglet_id)/2) ;ElasticLift_wing(1+length(ElasticLift_wing)/2:length(ElasticLift_wing))  ; ElasticLift_winglet(1+length(winglet_id)/2:end)];
ElasticMoment_updated = [ElasticMoment_wing(1:length(ElasticMoment_wing)/2) ; ElasticMoment_winglet(1:length(winglet_id)/2) ;ElasticMoment_wing(1+length(ElasticMoment_wing)/2:length(ElasticMoment_wing))  ; ElasticMoment_winglet(1+length(winglet_id)/2:end)];




