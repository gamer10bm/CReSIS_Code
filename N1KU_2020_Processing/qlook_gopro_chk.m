gid = 1:2;

figure(5)
clf

plot(qlook_data.Longitude,qlook_data.Latitude,'rs','MarkerFaceColor','r'); 
hold on
plot(gopro_data.lon(gid),gopro_data.lat(gid),'bo','MarkerFaceColor','b');
hold off

xlabel('Longitude')
ylabel('Latitude')