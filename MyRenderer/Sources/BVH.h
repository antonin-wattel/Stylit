// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#pragma once

#include <vector>
#include <memory>
#include <utility>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "BoundingBox.h"
#include "Ray.h"
#include "Mesh.h"
#include "Scene.h"

class BVH {
public:
    BVH(const std::shared_ptr<Scene> scene);

    BVH(const BVH& bvh);

    virtual ~BVH();

    BVH& operator= (const BVH& bvh);

    inline bool isLeaf() const { return (m_left == nullptr); }

    inline const BoundingBox& bbox() const { return m_bbox; }

    inline const BVH* left() const { return m_left; }

    inline const BVH* right() const { return m_right; }

    inline size_t meshIndex() const { return m_meshIndex; }

    inline size_t triangleIndex() const { return m_triangleIndex; }

    void intersect(const Ray& r, std::vector<std::pair<size_t,size_t>>& candidateMeshTrianglePairs) const;
    
private:
    BVH(const std::shared_ptr<Scene> scene, std::vector<std::pair<size_t,size_t> >& indexPairSet);
    BVH(const std::shared_ptr<Scene> scene, std::vector<std::pair<size_t, size_t> >& indexPairSet, size_t begin, size_t end);
    std::vector<std::pair<size_t, size_t> > makeIndexPairSet(const std::shared_ptr<Scene> scene);

    BoundingBox m_bbox;
    BVH* m_left;
    BVH* m_right;
    size_t m_meshIndex;
    size_t m_triangleIndex;
};
