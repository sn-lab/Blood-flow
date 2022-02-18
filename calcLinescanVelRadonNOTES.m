
    % if I_sign==1
    %     thetaMax = theta(n)-90;
    %     velocity = -1*(1/tand(thetaMax))*Tfactor*Xfactor;
    % elseif I_sign==0
    %     thetaMax = (theta(n)-90)*-1;
    %     velocity = (1/tand(-1*thetaMax))*Tfactor*Xfactor;
    % else
    %     if theta(n)<90
    %         if n-1<1
    %             n=2;
    %         end
    %         thetaMax = theta(n-1)-90;
    %     else
    %         thetaMax = theta(n)-90;
    %     end
    %     velocity = (1/tand(-1*thetaMax))*Tfactor*Xfactor;
    % end
    

    % maxVel = 60;
    % minVel = -60;
    % tol = 0.1;
    % thetas = [atand(0:tol:maxVel), atand(minVel:tol:0)+180];
    % % OR
    % maxSpeed = 60;
    % thetaHalf = atand(0:tol:maxSpeed);
    % thetas = [thetaHalf, -fliplr(thetaHalf)+180];