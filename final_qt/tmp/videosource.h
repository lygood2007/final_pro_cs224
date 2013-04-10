#ifndef VIDEOSOURCE_H
#define VIDEOSOURCE_H

#include <vlc/vlc.h>
#include <vlc/libvlc.h>
#include <QString>
#include <QMutex>

class VideoSource
{
public:
    VideoSource();
    virtual ~VideoSource();

    void playFile(QString file);

    void updateInterface();
    void changeVolume(int newVolume);
    void changePosition(int newPosition);

    void prerender();
    void postrender();

private:
    bool _isPlaying;
//    libvlc_exception_t _vlcexcep;
    libvlc_instance_t *_vlcinstance;
    libvlc_media_player_t *_mp;
    libvlc_media_t *_m;
};

struct ctx
{
    uchar*                  pixels;
    QMutex*                 mutex;
};


#endif // VIDEOSOURCE_H
