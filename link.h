#ifndef LINK_H
#define LINK_H

#include "macro.h"
#include "random.h"

typedef struct Link
{
    complex link;
}Link;

// general functions
void PrintLink(struct Link *link);
inline void CopyLink(const struct Link const * from, struct Link * destination) {destination->link = from->link;};
inline void SetToZero(struct Link * link){link->link = 0.0;};
void StartLink(struct Link * link);
void CorrectLink(struct Link * link);

// Link operations
inline void LinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * link3){link3->link = link1->link * link2->link;};
inline void LinkAdjLink(const struct Link const * link1, const struct Link const * link2, struct Link * link3){link3->link = link1->link * conj(link2->link);};
inline void AdjLinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * link3){link3->link = conj(link1->link) * link2->link;};
inline void AdjLinkAdjLink(const struct Link const * link1, const struct Link const * link2, struct Link * link3){link3->link = conj(link1->link) * conj(link2->link);};

inline void SumLinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple){staple->link += link1->link * link2->link;};
inline void SumLinkAdjLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple){staple->link += link1->link * conj(link2->link);};
inline void SumAdjLinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple){staple->link += conj(link1->link) * link2->link;};
inline void SumAdjLinkAdjLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple){staple->link += conj(link1->link) * conj(link2->link);};

inline void SubtractLinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple){staple->link -= link1->link * link2->link;};

#endif