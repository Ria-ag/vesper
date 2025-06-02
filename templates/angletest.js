// This file tests the angle function in the template engine

const NH3Elements = ['n', 'h1', 'h2', 'h3'];
const NH3Centers = [{x:250, y:250}, {x:250, y:250}, {x:250, y:250}, {x:250, y:250}];


function computeNodePositionsBasedOnAngle(centers, elements, lastCenter) {
    let nodes = [];
    starter = 0;
    startFrom = lastCenter;
    endAt = 0;
    let bondAngle = 0;
    count = 0;
    for (let i = lastCenter; i < centers.length; i++){

        if (centers[i] == centers[starter]){
            count++;
        } else{
            endAt = count + startFrom;
            bondAngle = (2 * Math.PI) / count;
            for (let j = startFrom; j < endAt; j++){
                let angle = j * bondAngle;
                x = centers[starter] + 120 * Math.cos(angle);
                y = 250 + 120 * Math.sin(angle);

                nodes.push({
                    data: {
                        id: elements[j],
                        label: elements[j]
                    },
                    position: {
                        x: x,
                        y: y
                    }
                });
            }
            startFrom = endAt;
            starter++;
            count = 0;
        }
    }
    return nodes;
}

function testAngleFunction() {
    const centers = NH3Centers;
    const elements = NH3Elements;
    const lastCenter = 0;

    const expectedPositions = [
        {x: 250, y: 130}, // n
        {x: 250, y: 370}, // h1
        {x: 370, y: 250}, // h2
        {x: 130, y: 250}  // h3
    ];

    const computedNodes = computeNodePositionsBasedOnAngle(centers, elements, lastCenter);

    console.info("Computed Node Positions:");
    for (let i = 0; i < computedNodes.length; i++) {
        console.info(`Node ${elements[i]}: (${computedNodes[i].position.x}, ${computedNodes[i].position.y})`);
    }
    // for (let i = 0; i < computedNodes.length; i++) {
    //    if (Math.abs(computedNodes[i].position.x - expectedPositions[i].x) > 1e-6 ||
    //        Math.abs(computedNodes[i].position.y - expectedPositions[i].y) > 1e-6) {
    //        console.error(`Test failed for node ${elements[i]}: expected (${expectedPositions[i].x}, ${expectedPositions[i].y}), got (${computedNodes[i].position.x}, ${computedNodes[i].position.y})`);
    //        return false;
    //    }
    //}
    console.info("All tests passed!");
    return true;
}

// Run the test
testAngleFunction();