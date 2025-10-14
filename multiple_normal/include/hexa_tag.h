#ifndef HEXA_TAG_H
#define HEXA_TAG_H

#include <functional>
#include <cstddef>

// Enum to define different types of tags.
enum EleSurfType {
    QUAD_SIDE_01,
    QUAD_SIDE_12,
    QUAD_SIDE_23,
    QUAD_SIDE_30,
    TRI_SIDE_01,
    TRI_SIDE_12,
    TRI_SIDE_20,
    TRI_TOP,
    TRI_BOTTOM,
    QUAD_TOP,
    QUAD_BOTTOM
};

// Struct representing a tag with a unit_id, layer, and type.
struct HexaTag {
    size_t unit_id;  // Unique ID for the unit.
    short layer;          // Layer index.
    EleSurfType type;     // Type associated with the tag.

    // Constructor with default type set to TRI_BOTTOM.
    HexaTag() : unit_id(0), layer(0), type(TRI_BOTTOM) {}
    HexaTag(size_t unit_id, short layer, EleSurfType type = TRI_BOTTOM)
        : unit_id(unit_id), layer(layer), type(type) {
    }
    HexaTag(int unit_id, short layer, EleSurfType type = TRI_BOTTOM)
        : unit_id(static_cast<size_t>(unit_id)), layer(layer), type(type) {
    }

    // Overload equality operator to compare HexaTags.
    bool operator==(const HexaTag& other) const {
        return unit_id == other.unit_id &&
            layer == other.layer &&
            type == other.type;
    }
};

// Overload inequality operator for convenience.
inline bool operator!=(const HexaTag& lhs, const HexaTag& rhs) {
    return !(lhs == rhs);
}

// Custom hash function for HexaTag to be used in unordered_map or unordered_set.
struct HexaTagHash {
    std::size_t operator()(const HexaTag& tag) const {
        // Compute hash values for each member.
        std::size_t h1 = std::hash<std::size_t>()(tag.unit_id);
        std::size_t h2 = std::hash<short>()(tag.layer);
        std::size_t h3 = std::hash<int>()(static_cast<int>(tag.type));
        // Combine the hash values using XOR and bit shifting.
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

#endif // HEXA_TAG_H
